# Step 7: validation of the panel method by anomaly injection, choice of K
# and of the penalty, comparison with the Hampel filter. About 35 minutes on
# 6 cores. Self-contained: the reference segmentation is refitted here.
# Run from the project root: Rscript scripts/05_calibrate.R
source("R/config.R")
source("R/load.R")
source("R/features.R")
source("R/checks.R")
source("R/segment.R")
source("R/panel.R")
source("R/classify.R")
source("R/hampel.R")
source("R/simulate.R")
suppressPackageStartupMessages(library(ggplot2))

for (d in c(config$figures_dir, config$tables_dir, config$results_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}
save_fig <- function(plot, name, width = 7, height = 4.5) {
  ggsave(file.path(config$figures_dir, paste0(name, ".png")), plot,
         width = width, height = height, dpi = 300, bg = "white")
}
theme_set(theme_minimal(base_size = 10))
cores <- max(1L, min(6L, parallel::detectCores() - 2L))

bills <- load_bills(config$data_path)
big <- big_sites(bills, config$big_site_kwh)
seg_args <- function(K, penalty) {
  list(cost = "biweight", penalty = penalty, plateau_min = config$plateau_min_months, K = K)
}

# ---- Ingredients of the simulated panels ------------------------------------
ref <- panel_fit(bills, month_effect(bills, big), big,
                 seg_args(config$calibration_reference$K, config$calibration_reference$penalty))
parts <- simulation_parts(ref$bills, ref$breaks, big, config$level_shift_min)
n_months <- length(parts[[1]]$month)
resid_acf <- apply(vapply(parts, function(p) stats::acf(p$resid, lag.max = 6, plot = FALSE)$acf[2:7],
                          numeric(6)), 1, stats::median)
# Spread of the residuals in units of the noise scale of the segmentation.
ref_sigma <- tapply(ref$bills$sigma, ref$bills$site, `[`, 1)
noise_ratio <- stats::median(vapply(parts, function(p) stats::sd(p$resid), numeric(1)) /
                               ref_sigma[names(parts)])
n_true_breaks <- sum(lengths(lapply(parts, `[[`, "break_index")))

settings <- expand.grid(K = config$calibration_K, beta = config$calibration_beta)
both <- function(type, sizes) data.frame(type = type, size = rep(sizes, each = 2), sign = c(1, -1))
cells <- rbind(data.frame(type = "none", size = 0, sign = 0),
               both("spike", config$calibration_spike_sizes),
               both("level_break", config$calibration_break_sizes),
               both("temporary", config$calibration_episode_sizes),
               both("pair", config$calibration_episode_sizes))
cells$id <- seq_len(nrow(cells))
expected_type <- c(spike = "point", level_break = "break", temporary = "temporary",
                   pair = "compensated_pair")

# One replication of one cell for the settings `which_settings`: panel and
# Hampel results. The seed depends on the replication only, so every cell and
# setting sees the same background panels and start months.
run_task <- function(r, cell, which_settings) {
  set.seed(config$calibration_seed + r)
  b <- simulate_panel(parts, config$calibration_block)
  u <- stats::runif(length(parts))
  span <- if (cell$type %in% c("temporary", "pair")) config$calibration_episode_months else 1L
  start <- setNames(anomaly_start(u, n_months, config$edge_months, span), names(parts))
  if (cell$type != "none") {
    b <- inject_anomaly(b, cell$type, cell$size, cell$sign, start, config$calibration_episode_months)
  }
  effect <- month_effect(b, character())
  hz <- hampel_bills(b, config$hampel_half_width, config$hampel_threshold)$z
  hampel <- data.frame(rep = r, cell = cell$id, threshold = config$hampel_thresholds,
                       hits = if (cell$type == "none") {
                         vapply(config$hampel_thresholds, function(h) sum(abs(hz) > h), numeric(1))
                       } else {
                         hampel_hits(hz, b, cell$type, start, config$calibration_episode_months,
                                     config$hampel_thresholds)
                       }, n = if (cell$type == "none") nrow(b) else length(start))
  panel <- do.call(rbind, lapply(which_settings, function(k) {
    fl <- detect_panel(b, effect, settings$K[k], settings$beta[k], config)
    if (cell$type == "none") {
      return(data.frame(setting = k, rep = r, cell = cell$id, n = nrow(b),
                        detected = false_alerts(fl, parts, config$calibration_break_window),
                        type_ok = NA, as_spike = NA, type = NA))
    }
    sc <- score_injection(fl, b, cell$type, cell$sign, start, config$calibration_episode_months,
                          config$calibration_break_window, config$level_shift_min)
    types <- table(factor(sc$type))
    data.frame(setting = k, rep = r, cell = cell$id, n = nrow(sc), detected = sum(sc$detected),
               type_ok = sum(sc$type == expected_type[[cell$type]], na.rm = TRUE),
               as_spike = sum(sc$as_spike, na.rm = TRUE),
               type = paste(names(types), types, sep = "=", collapse = ";"))
  }))
  list(panel = panel, hampel = hampel)
}
run_cells <- function(cell_ids, which_settings) {
  tasks <- expand.grid(r = seq_len(config$calibration_reps), cell = cell_ids)
  out <- parallel::mclapply(seq_len(nrow(tasks)), function(i) {
    run_task(tasks$r[i], cells[tasks$cell[i], ], which_settings)
  }, mc.cores = cores)
  failed <- vapply(out, inherits, logical(1), "try-error")
  if (any(failed)) stop("Calibration task failed: ", out[[which(failed)[1]]])
  list(panel = do.call(rbind, lapply(out, `[[`, "panel")),
       hampel = do.call(rbind, lapply(out, `[[`, "hampel")))
}

# ---- Phase 1: every setting on the background, spikes and breaks ------------
phase1_cells <- cells$id[cells$type %in% c("none", "spike", "level_break")]
timing <- system.time(res1 <- run_cells(phase1_cells, seq_len(nrow(settings))))[["elapsed"]]

rates <- function(panel) {
  agg <- stats::aggregate(cbind(n, detected, type_ok, as_spike) ~ setting + cell, panel, sum,
                          na.action = stats::na.pass)
  agg$rate <- agg$detected / agg$n
  agg$type_rate <- agg$type_ok / agg$detected
  agg$as_spike_rate <- agg$as_spike / agg$n
  cbind(agg, cells[agg$cell, c("type", "size", "sign")])
}
cal <- rates(res1$panel)
grid <- settings
grid$fa <- cal$rate[cal$type == "none"][match(seq_len(nrow(settings)), cal$setting[cal$type == "none"])]
mean_rate <- function(type) {
  x <- cal[cal$type == type, ]
  tapply(x$rate, x$setting, mean)[as.character(seq_len(nrow(settings)))]
}
grid$spike <- mean_rate("spike")
grid$level_break <- mean_rate("level_break")
grid$score <- (grid$spike + grid$level_break) / 2
grid$admissible <- grid$fa <= config$calibration_fa_max
chosen <- choose_setting(grid$fa, grid$score, config$calibration_fa_max, config$calibration_tie)
grid$chosen <- seq_len(nrow(grid)) == chosen

# ---- Phase 2: temporary gaps and pairs at the chosen setting ----------------
phase2_cells <- cells$id[cells$type %in% c("temporary", "pair")]
timing <- timing + system.time(res2 <- run_cells(phase2_cells, chosen))[["elapsed"]]
detection <- rbind(cal[cal$setting == chosen & cal$type != "none", ],
                   rates(res2$panel))

# Hampel on the same panels: false alert rate and detection by threshold.
hampel <- rbind(res1$hampel, res2$hampel)
hampel <- stats::aggregate(cbind(hits, n) ~ cell + threshold, hampel, sum)
hampel$rate <- hampel$hits / hampel$n
hampel <- cbind(hampel, cells[hampel$cell, c("type", "size", "sign")])
hampel_fa <- hampel[hampel$type == "none", c("threshold", "rate")]
matched <- min(c(Inf, hampel_fa$threshold[hampel_fa$rate <= grid$fa[chosen]]))
hampel_compared <- c(default = config$hampel_threshold, matched = matched)

# ---- Sensitivity: alerts on the real data for every setting ------------------
periods <- check_periods(bills)
real_effect <- month_effect(bills, big)
sensitivity <- do.call(rbind, lapply(seq_len(nrow(settings)), function(k) {
  pf <- panel_fit(bills, real_effect, big, seg_args(settings$K[k], settings$beta[k]))
  cfg <- config
  cfg$segment_K <- settings$K[k]
  fl <- classify_bills(pf$bills, pf$breaks, periods, big, cfg)
  fl <- fl[!fl$big, ]
  sus <- fl[fl$suspect, ]
  data.frame(setting = k, suspect_rows = nrow(sus),
             suspect_bills = nrow(unique(sus[!is.na(sus$episode) | sus$type == "irregular_period",
                                             c("site", "month")])),
             lasting_changes = sum(sus$type == "lasting_change"),
             large_breaks = sum(abs(pf$breaks$change) > config$level_shift_min &
                                  !pf$breaks$site %in% big),
             small_breaks = sum(abs(pf$breaks$change) <= config$level_shift_min &
                                  !pf$breaks$site %in% big),
             sites_suspect = length(unique(sus$site)))
}))
grid <- cbind(grid, sensitivity[-1])

# ---- Figure -----------------------------------------------------------------
type_labels <- c(spike = "Écart ponctuel", temporary = "Écart temporaire (3 mois)",
                 pair = "Paire compensée", level_break = "Rupture de niveau")
method_labels <- c(panel = "Panel (réglage retenu)",
                   hampel = sprintf("Hampel, seuil %g", hampel_compared[["default"]]))
fig <- rbind(
  data.frame(method = "panel", detection[c("type", "size", "sign", "rate")]),
  data.frame(method = "hampel",
             hampel[hampel$threshold == hampel_compared[["default"]] & hampel$type != "none",
                    c("type", "size", "sign", "rate")])
)
fig$type <- factor(type_labels[fig$type], levels = type_labels)
fig$method <- factor(method_labels[fig$method], levels = method_labels)
fig$sign <- factor(ifelse(fig$sign > 0, "hausse (paire : sous-facturation d'abord)",
                          "baisse (paire : surfacturation d'abord)"))
p <- ggplot(fig, aes(100 * size, 100 * rate, colour = method, linetype = sign)) +
  geom_line() +
  geom_point(size = 1) +
  scale_x_log10(breaks = c(1, 2, 5, 10, 20, 50, 100)) +
  scale_colour_manual(values = c("#08519c", "#fd8d3c"), name = NULL) +
  scale_linetype_manual(values = c("solid", "dashed"), name = NULL) +
  facet_wrap(~type, ncol = 2) +
  labs(x = "Amplitude injectée (%, échelle log)", y = "Part des anomalies signalées (%)",
       caption = "Rupture : signalée par Hampel quand l'un de ses deux premiers mois dépasse le seuil.") +
  guides(colour = guide_legend(ncol = 1), linetype = guide_legend(ncol = 1)) +
  theme(legend.position = "bottom", legend.box = "horizontal")
save_fig(p, "calibration_detection", height = 5.5)
# Conclusion: only K = 6 with the MBIC penalty keeps false alerts under 1 %
# (0.9 %); K = 3, the step 5 value, gives 5.1 %, and a higher penalty raises
# false alerts at every K (fewer plateaus, so slow drifts leave them). At
# K = 6, spikes of 10 % and 20 % are flagged 59 % and 88 % of the time, 5 %
# ones 16 %; breaks of 20 % are found 94 % of the time, 10 % ones 40 %.
# Breaks are always typed as breaks, spikes as points in 92 % of the
# detections; 3-month temporary gaps and pairs are flagged (pairs of 10 %:
# 94 %) but read as two close lasting changes (right type in 1 % of pairs).
# Hampel at threshold 3 has a similar false alert rate (1.0 %) and flags
# fewer spikes (50 % at 10 %, 78 % at 20 %) and fewer 3-month gaps (45 %
# against 83 % at 20 %). On the real data, suspect bills drop from 173 to 48.

# ---- Saved tables and numbers for the report ---------------------------------
write.csv(grid, file.path(config$tables_dir, "calibration_grid.csv"), row.names = FALSE)
write.csv(detection, file.path(config$tables_dir, "calibration_detection.csv"), row.names = FALSE)

results <- list(
  config = config[grep("^calibration_|^hampel_thresholds$|^level_shift_min$", names(config))],
  timing_seconds = timing,
  cores = cores,
  n_sites = length(parts),
  n_true_breaks = n_true_breaks,
  resid_acf = resid_acf,
  noise_ratio = noise_ratio,
  settings = settings,
  grid = grid,
  chosen = grid[chosen, ],
  calibration = cal,
  detection = detection,
  hampel = hampel,
  hampel_fa = hampel_fa,
  hampel_compared = hampel_compared,
  type_labels = type_labels
)
saveRDS(results, file.path(config$results_dir, "calibration.rds"))

# The segmentation settings of R/config.R must be those the rule chose.
config_beta <- if (identical(config$segment_penalty, "MBIC")) 3 else config$segment_penalty
if (config$segment_K != grid$K[chosen] || config_beta != grid$beta[chosen]) {
  stop(sprintf("R/config.R must set segment_K = %g and segment_penalty = %g (rule of step 7)",
               grid$K[chosen], grid$beta[chosen]))
}
cat("Calibration written to", config$output_dir, "in", round(timing), "s\n")
