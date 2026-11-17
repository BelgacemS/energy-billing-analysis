# Step 9 bis: the Box-Jenkins approach of the course done by hand on 0698BR,
# the only site forecast with a seasonal model, and the residual check of the
# per-day naive on the ordinary sites. Months before test_start only.
# Run from the project root: Rscript scripts/07b_box_jenkins.R
source("R/config.R")
source("R/load.R")
source("R/features.R")
source("R/models.R")
suppressPackageStartupMessages({
  library(ggplot2)
  library(forecast)
})

for (d in c(config$figures_dir, config$tables_dir, config$results_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

bills <- load_bills(config$data_path)
pre <- bills[bills$month < config$test_start, ]  # test months are never read here
site <- config$arima_season_sites
stopifnot(length(site) == 1)
d <- pre[pre$site == site, ]
y <- site_ts(d[order(d$month), ])  # log kWh per calendar day

# ---- Identification -----------------------------------------------------------
# Seasonal difference first (course 3), then the simple difference if the
# seasonally differenced series still wanders.
# Seasonal differencing tests on these months only (step 3 bis ran them on
# the 35 months, test months included).
nsdiffs_pre <- c(seas = nsdiffs(y, test = "seas"), ocsb = nsdiffs(y, test = "ocsb"))
d12 <- diff(y, lag = 12)
dd <- diff(d12)
acf_d12 <- acf(d12, plot = FALSE, lag.max = 3)$acf[2:4]
acf_dd <- acf(dd, plot = FALSE, lag.max = 12)$acf[-1]   # lags 1..12
pacf_dd <- pacf(dd, plot = FALSE, lag.max = 12)$acf
bound <- 1.96 / sqrt(length(dd))

png(file.path(config$figures_dir, "box_jenkins_0698BR.png"), width = 7, height = 4.5,
    units = "in", res = 300)
# ggtsdisplay draws its three panels itself (see 01b_course_view.R).
ggtsdisplay(dd, main = paste(site, "après une différence saisonnière et une différence simple"),
            lag.max = 12, theme = theme_minimal(base_size = 10))
invisible(dev.off())

# ---- Estimation: small orders around the identified model, same d and D ------
orders <- expand.grid(p = 0:1, q = 0:1, P = 0:1, Q = 0:1)
orders <- orders[orders$P + orders$Q <= 1, ]
candidates <- do.call(rbind, lapply(seq_len(nrow(orders)), function(k) {
  o <- orders[k, ]
  fit <- tryCatch(Arima(y, order = c(o$p, 1, o$q), seasonal = c(o$P, 1, o$Q)),
                  error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  lb <- arima_ljung_box(fit)
  data.frame(model = arima_label(fit), n_par = o$p + o$q + o$P + o$Q, aicc = fit$aicc,
             lb_lag = lb[["lag"]], lb_p = lb[["p_value"]], stringsAsFactors = FALSE)
}))
candidates <- candidates[order(candidates$aicc), ]
rownames(candidates) <- NULL

# ---- Validation of the retained model and comparison with auto.arima --------
best <- Arima(y, order = c(0, 1, 0), seasonal = c(0, 1, 0))
stopifnot(arima_label(best) == candidates$model[1])
lb_best <- arima_ljung_box(best)
auto_pre <- auto.arima(y)

# ---- Residuals of the per-day naive, ordinary sites ----------------------------
# The naive residuals are the monthly changes of log kWh per day; Ljung-Box
# with the lag rule of arima_ljung_box (no estimated parameter).
big <- big_sites(pre, config$big_site_kwh)
naive_lb <- do.call(rbind, lapply(split(pre[!pre$site %in% big, ], pre$site[!pre$site %in% big]),
                                  function(s) {
  r <- diff(site_ts(s[order(s$month), ]))
  h <- min(10, floor(length(r) / 5))
  bt <- Box.test(r, lag = h, type = "Ljung-Box")
  data.frame(site = s$site[1], lag = h, p_value = bt$p.value, acf1 = acf(r, plot = FALSE)$acf[2])
}))
naive_reject <- mean(naive_lb$p_value < config$ljung_box_alpha)

# Conclusion (29 months, January 2019 to May 2021): the seas test asks for a
# seasonal difference, the OCSB test does not (both do on 35 months). The seasonally differenced
# series still wanders (ACF 0.67, 0.51, 0.41 at lags 1 to 3, ndiffs 1). After
# both differences, 16 values are left: no ACF or PACF value leaves the
# +-0.49 bounds, so no AR or MA term is identified (with so few values, only
# a large autocorrelation would show). ARIMA(0,1,0)(0,1,0)[12] has the lowest
# AICc (-71.7, next -69.4) and white residuals (Ljung-Box on 5 lags,
# p = 0.93); auto.arima picks the same model. Its forecast is November plus
# the November to December change of the year before. Refitted on the 35
# months (09_forecast.R), auto.arima keeps the same order. Per-day naive on
# the 90 ordinary sites: Ljung-Box on 5 lags rejects at 5 % for 18 sites (one
# in five, 4.5 expected by chance), with a positive lag-1 autocorrelation
# (median 0.40): changes spread over two months. The level is close to a
# random walk, not exactly one.

write.csv(candidates, file.path(config$tables_dir, "box_jenkins_0698BR.csv"), row.names = FALSE)
saveRDS(list(
  site = site, n_months = length(y), nsdiffs_pre = nsdiffs_pre, n_dd = length(dd), bound = bound,
  acf_d12 = acf_d12, ndiffs_d12 = ndiffs(d12), acf_dd = acf_dd, pacf_dd = pacf_dd,
  candidates = candidates, retained = arima_label(best), lb_lag = lb_best[["lag"]],
  lb_p = lb_best[["p_value"]], sigma = sqrt(best$sigma2), auto_pre = arima_label(auto_pre),
  naive_lb = naive_lb, naive_reject = naive_reject,
  naive_acf1_median = median(naive_lb$acf1[naive_lb$p_value < config$ljung_box_alpha])
), file.path(config$results_dir, "box_jenkins.rds"))
cat("Box-Jenkins check written to", config$output_dir, "\n")
