# Step 3 bis: the course's first look (plot, season, ACF/PACF, STL,
# differencing, transformation) on representative series.
# Run from the project root: Rscript scripts/01b_course_view.R
source("R/config.R")
source("R/load.R")
source("R/features.R")
source("R/hampel.R")
suppressPackageStartupMessages({
  library(ggplot2)
  library(forecast)
})

for (d in c(config$figures_dir, config$tables_dir, config$results_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}
fig_path <- function(name) file.path(config$figures_dir, paste0(name, ".png"))
save_fig <- function(plot, name, width = 7, height = 4.5) {
  ggsave(fig_path(name), plot, width = width, height = height, dpi = 300, bg = "white")
}
theme_set(theme_minimal(base_size = 10))

bills <- load_bills(config$data_path)
big <- big_sites(bills, config$big_site_kwh)
small_bills <- bills[!bills$site %in% big, ]
small_mat <- log_kwh_matrix(small_bills)  # (90 sites, 35 months)
start <- as.integer(format(min(bills$month), c("%Y", "%m")))
as_ts <- function(x) ts(x, start = start, frequency = config$season_period)

# ---- Representative series (kWh per calendar day) --------------------------
typical <- pick_typical_sites(small_bills, small_mat, config$level_shift_window,
                              config$level_shift_min, config$hampel_half_width,
                              config$hampel_threshold)
site_series <- function(s) as_ts(bills$kwh_per_day[bills$site == s])
series <- list(
  total = as_ts(as.numeric(tapply(small_bills$kwh_per_day, small_bills$month, sum)))
)
for (s in big) series[[s]] <- site_series(s)
for (k in names(typical)) series[[typical[[k]]]] <- site_series(typical[[k]])
# "stable" is picked as the least noisy site, so it is labelled that way.
labels <- c(total = paste0("Total (", nrow(small_mat), " sites)"),
            setNames(paste(big, "(gros site)"), big),
            setNames(paste0(typical, " (", c("moins bruité", "rupture", "pic"), ")"), typical))

# ---- 1. Plot the series ----------------------------------------------------
mts <- do.call(cbind, series)
colnames(mts) <- labels[names(series)]
p <- autoplot(mts, facets = TRUE) +
  labs(x = NULL, y = "kWh par jour") +
  theme(strip.text = element_text(hjust = 0))
save_fig(p, "course_series", height = 8)
# Conclusion: the total of the ordinary sites rises by 7 % from 2019 to 2021
# with a summer bump each year; 0698BR combines growth and a strong summer
# peak; 0998BR is flat, then drops from May 2021. The typical sites show the
# three "odd" cases of the course: 0615BB, the least noisy site, still steps
# down by 8.6 % in February 2020, below the 10 % threshold of step 3 but about
# 30 times its month-to-month noise (0.27 %); 1240BR drops by 95 % in May 2021;
# 0274BR has one high bill in June 2021 and a low one in August.

# ---- 2. Season plots (total and big sites) ---------------------------------
for (s in c("total", big)) {
  p <- ggseasonplot(series[[s]], year.labels = TRUE, year.labels.left = TRUE,
                    season.labels = c("janv.", "févr.", "mars", "avr.", "mai", "juin",
                                      "juil.", "août", "sept.", "oct.", "nov.", "déc.")) +
    labs(title = NULL, x = NULL, y = "kWh par jour")
  save_fig(p, paste0("seasonplot_", s), width = 6, height = 4)
}
# Conclusion: only 0698BR repeats the same shape every year (trough in winter,
# peak in July-August), shifted upwards by its growth. The curves of the total
# are mostly shifted by the level; they share only a summer rise, whose peak
# moves between June and August.

# ---- 3. Variance: transformation or not ------------------------------------
# Within a series, as in the course: standard deviation against mean by year.
yearly <- do.call(rbind, lapply(c("total", big), function(s) {
  y <- series[[s]]
  yr <- floor(time(y))
  data.frame(series = s, year = sort(unique(yr)),
             mean = as.numeric(tapply(y, yr, mean)),
             sd = as.numeric(tapply(y, yr, stats::sd)))
}))
# Across sites: does the spread grow with the level? The noise is measured
# directly on kWh per day (MAD of first differences / sqrt(2)), so the test
# does not rely on the log it is meant to justify.
kpd_mat <- exp(small_mat)  # (90 sites, 35 months), kWh per day
site_stats <- data.frame(
  site = rownames(small_mat),
  mean = as.numeric(tapply(small_bills$kwh_per_day, small_bills$site, mean)[rownames(small_mat)]),
  sd = as.numeric(tapply(small_bills$kwh_per_day, small_bills$site, stats::sd)[rownames(small_mat)]),
  noise_abs = apply(kpd_mat, 1, function(x) stats::mad(diff(x)) / sqrt(2))
)
slope_sd <- sd_mean_slope(site_stats$mean, site_stats$sd)
slope_noise <- sd_mean_slope(site_stats$mean, site_stats$noise_abs)
p <- ggplot(site_stats, aes(mean, noise_abs)) +
  geom_point(colour = "steelblue") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "grey40", linewidth = 0.5) +
  scale_x_log10() + scale_y_log10() +
  labs(x = "kWh par jour, moyenne du site (échelle log)",
       y = "Bruit d'un mois à l'autre, en kWh par jour (échelle log)")
save_fig(p, "noise_vs_mean", width = 6, height = 4.5)
lambdas <- vapply(series, function(y) BoxCox.lambda(y), numeric(1))
# Conclusion: the yearly check of the course is inconclusive here. With three
# years, the mean of the total moves by 7 % only, and a yearly standard
# deviation also holds the trend and the season of that year; for 0698BR the
# standard deviation rises with the mean (coefficient of variation 0.08 to
# 0.11). The decision rests on the comparison across sites. The standard
# deviation of a whole series grows less than the level (slope 0.60, 95 % CI
# 0.27 to 0.93) because it also holds the level changes, which are not tied to
# the level; the month-to-month noise, which leaves them out, grows in
# proportion (slope 0.93, CI 0.66 to 1.21, which contains 1): the log puts all
# sites on the same relative scale. BoxCox.lambda (Guerrero) hits the bounds
# of its search interval (-1 or 2) on every series: on 35 months it compares
# only two yearly blocks, and as the course notes, the plot decides.

# ---- 4. ACF and PACF of the log series -------------------------------------
for (s in names(series)) {
  png(fig_path(paste0("tsdisplay_", s)), width = 7, height = 5, units = "in", res = 300)
  # ggtsdisplay draws its three panels itself; printing the result would add
  # a second page holding only the PACF.
  ggtsdisplay(log(series[[s]]), main = labels[[s]], lag.max = 2 * config$season_period,
              theme = theme_minimal(base_size = 10))
  dev.off()
}
# The course identifies a seasonal model on the differenced series.
png(fig_path("tsdisplay_total_diff"), width = 7, height = 5, units = "in", res = 300)
ggtsdisplay(diff(log(series$total)), main = paste(labels[["total"]], "après une différence"),
            lag.max = 2 * config$season_period, theme = theme_minimal(base_size = 10))
dev.off()
# Conclusion: on the log total the ACF decays slowly from 0.88 at lag 1 with
# no peak at lag 12 (0.02), and the PACF has a single spike at lag 1: the level
# wanders, as for a random walk. After one difference the lag-1
# autocorrelation is -0.07, far from the over-differencing sign of the course;
# the largest autocorrelation is then at lag 12 (0.29), but below the 95 %
# bound (0.34): a trace of the summer season, not significant on 34 points. For 1240BR the slow decay comes only
# from the break: a difference there would model the break, not a trend.

# ---- 5. STL decomposition, strength of trend and season --------------------
for (s in c("total", big)) {
  fit <- stl(log(series[[s]]), s.window = config$stl_s_window, robust = config$stl_robust)
  save_fig(autoplot(fit) + labs(title = labels[[s]], x = NULL), paste0("stl_", s), height = 5.5)
}
strength <- t(vapply(series, function(y) {
  stl_strength(log(y), config$stl_s_window, config$stl_robust)
}, numeric(2)))
big_season <- vapply(big, function(s) {
  fit <- stl(log(series[[s]]), s.window = config$stl_s_window, robust = config$stl_robust)
  exp(range(fit$time.series[, "seasonal"])) - 1
}, numeric(2))  # (min, max) x big sites, relative to the trend
# Same measure on every ordinary site, to place the representative series.
site_strength <- t(apply(small_mat, 1, function(x) {
  stl_strength(as_ts(x), config$stl_s_window, config$stl_robust)
}))
# Conclusion: with a fixed season, the STL season strength is 0.96 for 0698BR
# (seasonal component from -13 % to +10 % around the trend) and 0.72 for the
# total, above the 0.64 threshold of nsdiffs: the total has a regular season,
# small in amplitude (step 3: -1.5 % in winter to +3.2 % in August) but large
# next to its remainder, since adding 90 sites averages out their noise. Site
# by site the season is weaker: 67 of the 90 ordinary sites are below 0.64
# (median 0.46). The trend (level) strength has a median of 0.85: the level
# and its changes dominate.

# ---- 6. Differencing --------------------------------------------------------
# nsdiffs' default test measures the season with mstl, whose season may evolve
# from year to year (s.window = 11): on 35 months it also absorbs level
# changes. Its strength and the OCSB test are kept next to the fixed-season
# STL strength, and D = 1 is read as reliable only when they agree.
tests <- data.frame(
  series  = names(series),
  label   = unname(labels[names(series)]),
  ndiffs  = vapply(series, function(y) ndiffs(log(y)), numeric(1)),
  ndiffs_after_diff = vapply(series, function(y) ndiffs(diff(log(y))), numeric(1)),
  nsdiffs_seas = vapply(series, function(y) nsdiffs(log(y)), numeric(1)),
  nsdiffs_ocsb = vapply(series, function(y) nsdiffs(log(y), test = "ocsb"), numeric(1)),
  season_strength_mstl = vapply(series, function(y) mstl_season_strength(log(y)), numeric(1)),
  season_strength = strength[, "season"],
  trend_strength  = strength[, "trend"],
  lambda  = lambdas,
  row.names = NULL
)
# Conclusion: ndiffs asks for one difference on every series and none after
# it. A seasonal difference is supported only for 0698BR, where the seas test
# and the OCSB test agree; for the total the seas test (season allowed to
# evolve) says 1 and OCSB 0, and for 1240BR OCSB says 1, probably because of
# the break.
# This matches step 3: a weak shared season, a strong own season for 0698BR.
# For the total the evidence is mixed: both STL strengths are above the
# threshold and the seas test says D = 1, while OCSB says 0 and the lag-12
# autocorrelation of the differenced series stays below its bound. The season
# is real but small; rather than a seasonal difference, which would use 12 of
# the 35 months, it is handled as a month effect shared by the ordinary sites
# (steps 5 and 10).

# The same disagreement site by site: auto.arima picks D with the seas test by
# default, so it matters for step 9.
site_D <- data.frame(
  seas = apply(small_mat, 1, function(x) nsdiffs(as_ts(x))),
  ocsb = apply(small_mat, 1, function(x) nsdiffs(as_ts(x), test = "ocsb"))
)
site_D_table <- table(seas = site_D$seas, ocsb = site_D$ocsb)
# Conclusion: with its default test, auto.arima would difference seasonally
# 22 of the 90 ordinary sites and OCSB 10; they ask for it together on 6 and
# disagree on 20. Site by site, a seasonal difference would also copy each
# site's own December 2020 effect, which does not repeat (step 3). The tests do not settle it: step 9 compares auto.arima with its
# default and without season (D = 0) on the validation period.

# ---- Saved tables and numbers for the report -------------------------------
write.csv(tests, file.path(config$tables_dir, "course_view_tests.csv"), row.names = FALSE)
write.csv(yearly, file.path(config$tables_dir, "yearly_sd_mean.csv"), row.names = FALSE)
write.csv(cbind(site = rownames(site_strength), as.data.frame(site_strength)),
          file.path(config$tables_dir, "site_stl_strength.csv"), row.names = FALSE)
acf_at <- function(y, lags) as.numeric(stats::acf(y, lag.max = max(lags), plot = FALSE)$acf[lags + 1])
acf_total <- c(acf_at(log(series$total), c(1, 12)), acf_at(diff(log(series$total)), c(1, 12)),
               1.96 / sqrt(length(series$total) - 1))
names(acf_total) <- c("lag1", "lag12", "lag1_after_diff", "lag12_after_diff", "bound_after_diff")
stable_shift <- max_level_shift(small_mat[typical[["stable"]], , drop = FALSE], config$level_shift_window)
stable_noise <- site_noise(small_mat[typical[["stable"]], , drop = FALSE])

# Same calendar months in the first and last year (the last one is partial).
tot <- series$total
last_months <- cycle(tot)[floor(time(tot)) == max(floor(time(tot)))]
in_year <- function(y) floor(time(tot)) == y & cycle(tot) %in% last_months
growth_total <- mean(tot[in_year(max(floor(time(tot))))]) / mean(tot[in_year(min(floor(time(tot))))]) - 1

results <- list(
  typical          = typical,
  site_D_table     = site_D_table,
  acf_total        = acf_total,
  stable_shift     = exp(stable_shift$shift) - 1,
  stable_shift_month = stable_shift$shift_month,
  stable_noise     = exp(stable_noise) - 1,
  big_season       = big_season,
  growth_total     = growth_total,
  labels           = labels,
  tests            = tests,
  yearly           = yearly,
  slope_sd         = slope_sd,
  slope_noise      = slope_noise,
  site_strength_q  = apply(site_strength, 2, stats::quantile, probs = c(0.1, 0.5, 0.9)),
  n_season_weak    = sum(site_strength[, "season"] < config$nsdiffs_seas_threshold)
)
saveRDS(results, file.path(config$results_dir, "course_view.rds"))
cat("Course view written to", config$output_dir, "\n")
