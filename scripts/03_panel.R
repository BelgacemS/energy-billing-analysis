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
  list(cost = cost, penalty = config$segment_penalty, minseglen = config$segment_minseglen,
       plateau_min = config$plateau_min_months, nquantiles = config$segment_nquantiles)
}

# ---- Synthetic check: does each cost separate spikes from breaks? ----------
synthetic <- do.call(rbind, lapply(c("normal", "np"), function(cost) {
  transform(synthetic_segmentation(
    n = length(unique(bills$month)), reps = config$synthetic_reps,
    seed = config$synthetic_seed, step_sizes = config$synthetic_step_sizes,
    spike_sizes = config$synthetic_spike_sizes, step_month = config$synthetic_step_month,
    spike_month = config$synthetic_spike_month, seg_args = seg_args(cost)), cost = cost)
}))
synthetic <- synthetic[order(synthetic$cost, synthetic$noise, synthetic$scenario,
                             synthetic$size), ]
# Conclusion: the Normal cost with one-month segments separates spikes from
# breaks. On a constant level it finds a lasting break in 5.5 % of the runs
# (Gaussian noise) and 6.5 % (t3); a spike of 3 to 8 noise units does not raise
# that rate (3 to 5.5 %), and a spike of 8 units is set apart as an excursion
# in 94 % of the runs. A step of 3 units is always found, at the exact month in
# 78 % (Gaussian) to 88 % (t3) of the runs; a step of 1 unit in 40 to 63 %.
# The nonparametric cost (cpt.np, same MBIC penalty) never sets a spike apart
# (at most 5 %) and finds a false break on a constant level in 28 to 36 % of
# the runs: an empirical distribution cannot be estimated on one month. The
# Normal cost is kept, no fallback is needed; the penalty is set at step 7.

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
shifts$has_break <- shifts$site %in% breaks$site
shifts$has_big_break <- shifts$site %in% breaks$site[abs(breaks$change) > config$level_shift_min]
agreement <- table(level_change = shifts$level_change, has_break = shifts$has_break)

# ---- Figures ---------------------------------------------------------------
# Sites cited in the subject, plus 1240BR (largest drop).
example_sites <- c("1240BR", "0089KO", "1067KO", "1088KO")
ex <- pb[pb$site %in% example_sites, ]
ex$site <- factor(ex$site, levels = example_sites)
p <- ggplot(ex, aes(month)) +
  geom_step(aes(y = exp(level)), colour = "grey55", direction = "mid") +
  geom_line(aes(y = exp(expected_log)), colour = "steelblue") +
  geom_point(aes(y = exp(x), colour = excursion), size = 1) +
  scale_colour_manual(values = c(`FALSE` = "black", `TRUE` = "firebrick"),
                      labels = c("facture", "excursion courte"), name = NULL) +
  scale_y_log10(labels = scales::label_number(big.mark = " ")) +
  scale_x_date(date_breaks = "6 months", date_labels = "%m/%Y") +
  facet_wrap(~site, ncol = 2, scales = "free_y") +
  labs(x = NULL, y = "kWh par jour (échelle log)",
       caption = "Gris : palier. Bleu : valeur attendue (palier + effet du mois).") +
  theme(legend.position = "bottom")
save_fig(p, "panel_examples", height = 5)
# Conclusion: the expected value follows the plateaus and keeps the isolated
# bills apart. 1240BR drops 95 % in May 2021, a lasting break (its noisier last
# months form one plateau). 0089KO bills 201 kWh in June 2019 and 5,133 in
# August against about 1,500 (z = -43 and +25), then steps up 59 % in
# September. 1067KO drops 54 % in July 2019 after two intermediate months, and
# 1088KO triples in August 2020 after a two-month ramp: both ramps are short
# excursions dated before the new plateau.

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
# Conclusion: with the default MBIC penalty, all 90 ordinary sites have at least
# one lasting break (311 in all, 225 of 10 % or less), 61 sites one above 10 %.
# Of the 67 sites with a level change above 10 % at step 3, 58 have a break
# above 10 % here; 3 of the 23 others do. The small breaks come from slow
# drifts and site-specific summers: the noise scale comes from first
# differences and the plateau residuals are autocorrelated (median lag-1
# 0.24), so the penalty is too lenient; step 7 sets it. Breaks cluster in June
# 2019 (26 sites), April 2020 (24), August 2021 (22) and December 2019 (20);
# above 10 %, April 2020 (10), December 2019 (9) and June 2019 (8).

# ---- Saved tables and numbers for the report -------------------------------
write.csv(synthetic, file.path(config$tables_dir, "segment_synthetic.csv"), row.names = FALSE)
write.csv(breaks, file.path(config$tables_dir, "panel_breaks.csv"), row.names = FALSE)
write.csv(effect, file.path(config$tables_dir, "panel_month_effect.csv"), row.names = FALSE)

small_sites <- setdiff(unique(bills$site), big)
results <- list(
  config = config[c("segment_cost", "segment_penalty", "segment_minseglen",
                    "plateau_min_months", "synthetic_reps", "synthetic_step_sizes",
                    "synthetic_spike_sizes")],
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
  shifts          = shifts,
  n_excursions    = sum(pb$excursion),
  n_excursion_sites = length(unique(pb$site[pb$excursion])),
  sigma_range     = range(tapply(pb$sigma, pb$site, unique)[small_sites]),
  sigma_median    = median(tapply(pb$sigma, pb$site, unique)[small_sites]),
  example_sites   = example_sites
)
saveRDS(results, file.path(config$results_dir, "panel.rds"))
cat("Panel written to", config$output_dir, "\n")
