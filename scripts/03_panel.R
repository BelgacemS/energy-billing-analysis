# Step 5: expected value of each bill by a robust panel, and level breaks.
# Run from the project root: Rscript scripts/03_panel.R
source("R/config.R")
source("R/load.R")
source("R/features.R")
source("R/segment.R")
source("R/panel.R")
suppressPackageStartupMessages(library(ggplot2))

for (d in c(config$figures_dir, config$tables_dir, config$results_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}
save_fig <- function(plot, name, width = 7, height = 4.5) {
  ggsave(file.path(config$figures_dir, paste0(name, ".png")), plot,
         width = width, height = height, dpi = 300, bg = "white")
}
theme_set(theme_minimal(base_size = 10))

bills <- load_bills(config$data_path)
big <- big_sites(bills, config$big_site_kwh)
seg_args <- function(cost) {
  list(cost = cost, penalty = config$segment_penalty, plateau_min = config$plateau_min_months,
       K = config$segment_K, nquantiles = config$segment_nquantiles)
}

# ---- Synthetic check: which cost separates anomalies from breaks? ---------
# The three costs share the same minimal plateau, threshold K and MBIC
# penalty; only the loss differs.
synthetic <- do.call(rbind, lapply(c("biweight", "normal", "np"), function(cost) {
  transform(synthetic_segmentation(
    n = length(unique(bills$month)), reps = config$synthetic_reps,
    seed = config$synthetic_seed, step_sizes = config$synthetic_step_sizes,
    spike_sizes = config$synthetic_spike_sizes, step_month = config$synthetic_step_month,
    anomaly_month = config$synthetic_anomaly_month, seg_args = seg_args(cost)), cost = cost)
}))
synthetic <- synthetic[order(synthetic$cost, synthetic$noise, synthetic$scenario,
                             synthetic$size), ]
# Conclusion: with the same minimal plateau, threshold and penalty, only the
# biweight cost separates anomalies from breaks. Its false break rate stays at
# 3 to 7.5 % (Gaussian noise) whether the level is constant or has a spike or
# a compensated pair, and at 2 to 9.5 % with t3 noise; a spike or a pair of 8
# noise units is flagged in 94 to 100 % of the runs. The Gaussian cost turns a
# pair of 5 to 8 units into a false plateau in 90 to 100 % of the runs, and a
# spike of 8 units in 46 % (t3: 80 %). The nonparametric cost finds a false
# break in 26 to 36 % of the runs whatever the scenario. All three find a step
# of 3 units, at the exact month in 78 to 88 % of the runs. Biweight kept.

# ---- Panel on the data -----------------------------------------------------
# Month effect from the trend of the classical decomposition, estimated apart
# from the segmentation. Re-estimating it on the plateaus (backfitting) was
# tried and dropped: the lenient penalty turns summers into plateaus that
# absorb the season (summer effect shrunk, no convergence after one step).
effect <- month_effect(bills, big)
panel <- panel_fit(bills, effect, big, seg_args(config$segment_cost))
pb <- panel$bills
breaks <- panel$breaks
breaks$big <- breaks$site %in% big

# Lasting breaks against the level change measure of step 3 (ordinary sites).
mat <- log_kwh_matrix(bills)
shifts <- max_level_shift(mat[!rownames(mat) %in% big, ], config$level_shift_window)
shifts$level_change <- abs(shifts$shift) > log(1 + config$level_shift_min)
shifts$has_big_break <- shifts$site %in% breaks$site[abs(breaks$change) > config$level_shift_min]
agreement <- table(level_change = shifts$level_change, big_break = shifts$has_big_break)

# Criterion of step 5: the 1240BR drop of May 2021 is a lasting break.
b1240 <- breaks[breaks$site == "1240BR" & breaks$month == as.Date("2021-05-01"), ]
stopifnot(nrow(b1240) == 1, b1240$change < -0.5)

# Lag-1 autocorrelation of the residuals around the plateaus (ordinary sites,
# excursions left out): above 0, the noise scale from first differences
# understates the spread around a plateau.
small_pb <- pb[!pb$site %in% big & !pb$excursion, ]
resid_acf1 <- vapply(split(small_pb$resid, small_pb$site),
                     function(r) stats::acf(r, lag.max = 1, plot = FALSE)$acf[2], numeric(1))

# ---- Figures ---------------------------------------------------------------
# Sites cited in the subject, 1240BR (largest drop) and 1911BB (ramp after
# normal months, the case that ruled out the Gaussian two-pass version).
example_sites <- c("1240BR", "0089KO", "1911BB", "1088KO")
ex <- pb[pb$site %in% example_sites, ]
ex$site <- factor(ex$site, levels = example_sites)
p <- ggplot(ex, aes(month)) +
  geom_step(aes(y = exp(level)), colour = "grey55", direction = "mid") +
  geom_line(aes(y = exp(expected_log)), colour = "steelblue") +
  geom_point(aes(y = exp(x), colour = excursion), size = 1) +
  scale_colour_manual(values = c(`FALSE` = "black", `TRUE` = "firebrick"),
                      labels = c("facture", sprintf("hors palier (écart > %g fois le bruit)",
                                                     config$segment_K)), name = NULL) +
  scale_y_log10(labels = scales::label_number(big.mark = " ")) +
  scale_x_date(date_breaks = "6 months", date_labels = "%m/%Y") +
  facet_wrap(~site, ncol = 2, scales = "free_y") +
  labs(x = NULL, y = "kWh par jour (échelle log)",
       caption = "Gris : palier. Bleu : valeur attendue (palier + effet du mois).") +
  theme(legend.position = "bottom")
save_fig(p, "panel_examples", height = 5)
# Conclusion: the expected value follows the plateaus and leaves the anomalous
# bills apart. 1240BR drops 94 % in May 2021, a lasting break; its noisier
# last months form two plateaus. 0089KO bills 201 kWh in June 2019 (z = -43)
# and 5,133 in August (z = +15), then stays on a level 60 % higher; the break
# is dated in August, one month early, as the catch-up and the new level
# cannot be told apart. 1911BB keeps its level until August 2019 and triples
# from October after a ramp in September (z = +91 and -40 for the ramp
# months). 1088KO triples from June 2020, its two ramp months lying below the
# new plateau.

br_small <- breaks[!breaks$big, ]
br_small$size <- ifelse(abs(br_small$change) > config$level_shift_min,
                        sprintf("plus de %d %%", round(100 * config$level_shift_min)),
                        sprintf("%d %% ou moins", round(100 * config$level_shift_min)))
p <- ggplot(br_small, aes(month, fill = size)) +
  geom_bar() +
  scale_fill_manual(values = c("grey70", "steelblue"), name = "Changement de niveau") +
  scale_x_date(date_breaks = "3 months", date_labels = "%m/%Y") +
  labs(x = "Premier mois du nouveau palier", y = "Nombre de sites ordinaires") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")
save_fig(p, "break_months")
# Conclusion: with the default MBIC penalty, all 90 ordinary sites have at
# least one lasting break (322 in all, 242 of 10 % or less), 58 sites one
# above 10 %. Of the 67 sites with a level change above 10 % at step 3, 56
# have a break above 10 % here; 2 of the 23 others do. The small breaks follow
# slow drifts and site-specific summers: the plateau residuals are
# autocorrelated (median lag-1 0.15, excursions left out), so the noise scale
# from first differences is too small and the penalty too lenient; step 7 sets
# K and the penalty. Breaks cluster in April 2020 (24 sites), June 2019 and
# August 2021 (21 each); above 10 %, April 2020 (9) and November 2019 (7).

# ---- Saved tables and numbers for the report -------------------------------
write.csv(synthetic, file.path(config$tables_dir, "segment_synthetic.csv"), row.names = FALSE)
write.csv(breaks, file.path(config$tables_dir, "panel_breaks.csv"), row.names = FALSE)
write.csv(effect, file.path(config$tables_dir, "panel_month_effect.csv"), row.names = FALSE)

small_sites <- setdiff(unique(bills$site), big)
results <- list(
  config = config[c("segment_cost", "segment_penalty", "segment_K", "plateau_min_months",
                    "synthetic_reps", "synthetic_step_sizes", "synthetic_spike_sizes",
                    "synthetic_step_month")],
  n_months        = length(unique(bills$month)),
  synthetic       = synthetic,
  effect          = effect,
  panel           = pb,
  breaks          = breaks,
  n_breaks_small  = nrow(br_small),
  n_sites_break   = length(unique(br_small$site)),
  n_sites_big_break = length(unique(br_small$site[abs(br_small$change) > config$level_shift_min])),
  n_small         = length(small_sites),
  break_months    = sort(table(format(br_small$month, "%Y-%m")), decreasing = TRUE),
  agreement       = agreement,
  resid_acf1_median = median(resid_acf1),
  shifts          = shifts,
  n_excursions    = sum(pb$excursion),
  n_excursion_sites = length(unique(pb$site[pb$excursion])),
  sigma_range     = range(tapply(pb$sigma, pb$site, unique)[small_sites]),
  sigma_median    = median(tapply(pb$sigma, pb$site, unique)[small_sites]),
  example_sites   = example_sites
)
saveRDS(results, file.path(config$results_dir, "panel.rds"))
cat("Panel written to", config$output_dir, "\n")
