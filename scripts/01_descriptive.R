# Step 3: descriptive statistics of the billing panel.
# Run from the project root: Rscript scripts/01_descriptive.R
source("R/config.R")
source("R/load.R")
source("R/features.R")
suppressPackageStartupMessages(library(ggplot2))

for (d in c(config$figures_dir, config$tables_dir, config$results_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}
save_fig <- function(plot, name, width = 7, height = 4.5) {
  ggsave(file.path(config$figures_dir, paste0(name, ".png")), plot,
         width = width, height = height, dpi = 300, bg = "white")
}
theme_set(theme_minimal(base_size = 10))
month_labels <- c("janv.", "févr.", "mars", "avr.", "mai", "juin",
                  "juil.", "août", "sept.", "oct.", "nov.", "déc.")

bills <- load_bills(config$data_path)
mat <- log_kwh_matrix(bills)  # (sites, months)
big <- big_sites(bills, config$big_site_kwh)
small_bills <- bills[!bills$site %in% big, ]

# ---- Levels and weight of the big sites ------------------------------------
sites <- data.frame(
  site       = rownames(mat),
  suffix     = substr(rownames(mat), 5, 6),
  median_kwh = as.numeric(tapply(bills$kwh, bills$site, median)[rownames(mat)]),
  total_kwh  = as.numeric(tapply(bills$kwh, bills$site, sum)[rownames(mat)]),
  stringsAsFactors = FALSE
)
sites$share <- sites$total_kwh / sum(sites$total_kwh)
sites$big <- sites$site %in% big
# Month-to-month volatility: MAD of the first differences of log kWh per day.
sites$volatility <- apply(mat, 1, function(x) stats::mad(diff(x)))

p <- ggplot(sites, aes(median_kwh, fill = big)) +
  geom_histogram(bins = 40, colour = "white") +
  scale_x_log10(labels = scales::label_number(big.mark = " ")) +
  scale_fill_manual(values = c(`FALSE` = "grey55", `TRUE` = "firebrick"),
                    labels = c("autres sites", "2 gros sites"), name = NULL) +
  labs(x = "kWh facturés par mois (médiane du site, échelle log)", y = "Nombre de sites")
save_fig(p, "levels")
# Conclusion: the median site bills 2,257 kWh a month; the 90 ordinary sites
# range from 544 to 6,236 kWh; 0698BR and 0998BR (433,000 and 347,000 kWh)
# carry 81 % of the kWh, so any total in kWh is driven by them and they are
# treated apart.

# ---- Month profile: raw kWh against kWh per calendar day -------------------
# Profiles use complete calendar years only: with the site-year median, the
# 11 months of 2021 put an exact zero on one month per site.
full_years <- complete_years(bills)
dev_day_all <- seasonal_deviation(small_bills, "kwh_per_day")
dev_day <- dev_day_all[dev_day_all$year %in% full_years, ]
dev_raw <- seasonal_deviation(small_bills, "kwh")
dev_raw <- dev_raw[dev_raw$year %in% full_years, ]
profile <- rbind(
  transform(common_profile(dev_raw), unit = "kWh du mois"),
  transform(common_profile(dev_day), unit = "kWh par jour")
)
p <- ggplot(profile, aes(moy, ratio, colour = unit)) +
  geom_hline(yintercept = 1, colour = "grey70") +
  geom_line() + geom_point() +
  scale_x_continuous(breaks = 1:12, labels = month_labels) +
  scale_y_continuous(breaks = seq(0.9, 1.06, by = 0.02),
                     labels = scales::label_percent(accuracy = 1)) +
  scale_colour_manual(values = c("grey55", "steelblue"), name = NULL) +
  labs(x = NULL, y = "Rapport à la médiane annuelle du site (médiane des sites)",
       caption = "Sites ordinaires, années complètes 2019 et 2020.")
save_fig(p, "month_profile")
# Conclusion: in raw kWh February is 9 % below the site's annual median and
# the other months follow their length (April and November -2.7 % and -1.6 %,
# July and August +3.7 % and +4.8 %). Per calendar day February is level with
# January and March (-1.5 to -1.7 %): its dip is a calendar effect. What
# remains is a small real seasonality, from -1.5 to -1.7 % in winter to +2.2
# to +3.2 % in July and August. Limit: the site-year median also absorbs part
# of the growth within a year, which pushes early months down and late months
# up.

# ---- Dispersion of the site-specific month effects -------------------------
eff <- specific_effect(dev_day)
site_eff <- stats::aggregate(effect ~ site + moy, eff, mean)
eff_lim <- 0.25
n_hidden <- sum(abs(exp(site_eff$effect) - 1) > eff_lim)
p <- ggplot(site_eff, aes(factor(moy), exp(effect) - 1)) +
  geom_hline(yintercept = 0, colour = "grey70") +
  geom_boxplot(outlier.size = 0.6, fill = "grey90") +
  scale_x_discrete(labels = month_labels) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
  coord_cartesian(ylim = c(-eff_lim, eff_lim)) +
  labs(x = NULL, y = "Effet propre du mois (écart au profil commun)",
       caption = sprintf("Moyenne 2019-2020 par site. %d valeurs sur %d hors de ±%d %% non affichées.",
                         n_hidden, nrow(site_eff), round(100 * eff_lim)))
save_fig(p, "month_effect_dispersion")
# Conclusion: around the common profile, 90 % of the site-specific month
# effects lie between -12 % and +14 %, about five times the amplitude of the
# common profile itself (-1.7 % to +3.2 %). Part of the tails comes from level
# changes within a year, not from seasonality (see the reproducibility check).

# ---- Reproducibility of the site-specific month effect ---------------------
# Robustness: same measure with a rolling 13-month level, less sensitive to
# level changes within a year.
dev_roll <- seasonal_deviation(small_bills, "kwh_per_day", centre = "rolling")
# The 2020-2021 pair keeps the 11 months of 2021 (see seasonal_deviation).
repro <- rbind(
  transform(effect_reproducibility(dev_day_all, 2019, 2020), centre = "year"),
  transform(effect_reproducibility(dev_day_all, 2020, 2021), centre = "year"),
  transform(effect_reproducibility(dev_roll, 2019, 2020), centre = "rolling"),
  transform(effect_reproducibility(dev_roll, 2020, 2021), centre = "rolling")
)
# Share of positive correlations, against 0 (Wilcoxon signed rank).
repro_test <- do.call(rbind, lapply(split(repro, repro$centre), function(r) {
  data.frame(centre = r$centre[1], n_pairs = nrow(r), n_positive = sum(r$rho > 0),
             median_rho = median(r$rho), wilcoxon_p = stats::wilcox.test(r$rho)$p.value)
}))
dec_rho_test <- function(d) {
  dd <- merge(d[d$year == 2019 & d$moy == 12, c("site", "dev")],
              d[d$year == 2020 & d$moy == 12, c("site", "dev")], by = "site")
  suppressWarnings(stats::cor.test(dd$dev.x, dd$dev.y, method = "spearman"))$p.value
}
dec_p <- c(year = dec_rho_test(dev_day_all), rolling = dec_rho_test(dev_roll))

# Points beyond the frame are drawn on its edge as triangles, so none is lost.
dec <- merge(eff[eff$year == 2019 & eff$moy == 12, c("site", "effect")],
             eff[eff$year == 2020 & eff$moy == 12, c("site", "effect")],
             by = "site", suffixes = c("_2019", "_2020"))
dec_lim <- 0.2
dec$x <- exp(dec$effect_2019) - 1
dec$y <- exp(dec$effect_2020) - 1
dec$edge <- abs(dec$x) > dec_lim | abs(dec$y) > dec_lim
dec$x <- pmax(pmin(dec$x, dec_lim), -dec_lim)
dec$y <- pmax(pmin(dec$y, dec_lim), -dec_lim)
p <- ggplot(dec, aes(x, y, shape = edge)) +
  geom_hline(yintercept = 0, colour = "grey70") +
  geom_vline(xintercept = 0, colour = "grey70") +
  geom_point(colour = "steelblue", alpha = 0.8) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), guide = "none") +
  scale_x_continuous(labels = scales::label_percent(accuracy = 1)) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
  labs(x = "Effet propre de décembre 2019", y = "Effet propre de décembre 2020",
       caption = sprintf("%d sites sur %d au-delà de ±%d %%, ramenés au bord (triangles).",
                         sum(dec$edge), nrow(dec), round(100 * dec_lim)))
save_fig(p, "december_reproducibility", width = 5, height = 4.5)
# Conclusion: a site's own December effect barely carries over from 2019 to
# 2020: Spearman rho -0.09 with the site-year centre, 0.24 (p = 0.02) with the
# rolling centre. Over the 23 month pairs (2019-2020 and 2020-2021), rho stays
# between -0.32 and +0.40, with a median of 0.06 (year) and 0.12 (rolling);
# with the rolling centre 17 of 23 are positive (Wilcoxon p = 0.01). There is
# a weak persistence, of the order of 0.1: a site's own month effect estimated
# on two years would be mostly noise, so the month effect is shared across
# sites (combined seasonal indices, Withycombe 1989).

# ---- Growth 2019 to 2021 and level changes ---------------------------------
growth <- site_growth(bills, 2019, 2021)
shifts <- max_level_shift(mat, config$level_shift_window)
sites <- merge(merge(sites, growth, by = "site"), shifts, by = "site")
sites$level_change <- abs(sites$shift) > log(1 + config$level_shift_min)

ord <- sites[order(sites$growth), ]
# Log scale on 1 + growth: a halving and a doubling are symmetric, and the
# extreme sites (-93 % to +218 %) stay visible.
growth_breaks <- c(-0.9, -0.5, -0.2, 0, 0.25, 1, 2)
p <- ggplot(ord, aes(seq_along(site), 1 + growth, colour = suffix)) +
  geom_hline(yintercept = 1, colour = "grey70") +
  geom_point() +
  scale_y_log10(breaks = 1 + growth_breaks,
                labels = scales::label_percent(accuracy = 1, style_positive = "plus")(growth_breaks)) +
  labs(x = "Sites, du plus en baisse au plus en hausse",
       y = "Évolution 2019 à 2021 (kWh par jour, janv. à nov.)", colour = "Suffixe")
save_fig(p, "growth")
# Conclusion: median growth of the ordinary sites is +7 % from 2019 to 2021,
# but the spread is wide (-93 % for 1240BR to +218 % for 1088KO). The extremes
# are jumps, not trends: 1088KO, 1978KO, 1986KO and 1093KO all step up in June
# 2020.

# ---- Heat map: deviation from the site level -------------------------------
site_order <- sites$site[order(sites$suffix, sites$growth)]
heat <- data.frame(
  site  = factor(rep(rownames(mat), ncol(mat)), levels = rev(site_order)),
  month = as.Date(paste0(rep(colnames(mat), each = nrow(mat)), "-01")),
  dev   = as.vector(mat - apply(mat, 1, median))
)
heat$dev_pct <- pmax(pmin(exp(heat$dev) - 1, 0.3), -0.3)
p <- ggplot(heat, aes(month, site, fill = dev_pct)) +
  geom_tile() +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick",
                       limits = c(-0.3, 0.3), labels = scales::label_percent(),
                       name = "Écart au\nniveau médian\n(borné à 30 %)") +
  scale_x_date(date_breaks = "6 months", date_labels = "%m/%Y", expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme(axis.text.y = element_text(size = 4), panel.grid = element_blank())
save_fig(p, "heatmap", width = 7, height = 9)
# Conclusion: most rows change colour in steps or fast drifts rather than
# slowly: 67 of the 90 ordinary sites (74 %) have a level change above 10 %
# between two 6-month windows (the measure does not separate a step from a
# drift of about 2 % a month; step 5 does). The change months cluster (July
# and October 2019, May-June 2020), which suggests network-wide events; July
# 2019 is also the first month the 6-month window can test, so it collects
# earlier changes. Cases cited in the
# subject show up: 1067KO (July 2019), 0254BB (end of 2019), 1088KO, 1978KO and
# 1986KO (June 2020), plus 1037KO down from November 2020 and 1240BR down 95 %
# from May 2021.

# ---- Suffixes --------------------------------------------------------------
small_sites <- sites[!sites$big, ]
suffix_tab <- do.call(rbind, lapply(split(small_sites, small_sites$suffix), function(s) {
  data.frame(suffix = s$suffix[1], n_sites = nrow(s),
             median_kwh = median(s$median_kwh),
             median_growth = median(s$growth),
             share_level_change = mean(s$level_change),
             median_volatility = median(s$volatility))
}))
# With about 30 sites per suffix, differences are tested before being read.
suffix_tests <- c(
  level_change_fisher_p = stats::fisher.test(table(small_sites$suffix, small_sites$level_change))$p.value,
  growth_kruskal_p = stats::kruskal.test(growth ~ factor(suffix), small_sites)$p.value
)
suffix_profile <- common_profile(dev_day, by = "suffix")
p <- ggplot(suffix_profile, aes(moy, ratio, colour = suffix)) +
  geom_hline(yintercept = 1, colour = "grey70") +
  geom_line() + geom_point() +
  scale_x_continuous(breaks = 1:12, labels = month_labels) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
  labs(x = NULL, y = "kWh par jour, rapport à la médiane annuelle", colour = "Suffixe")
save_fig(p, "suffix_profile")
# Conclusion: the per-day profiles of the three suffixes stay within about 2
# points of each other in every month (August peak +3.1 to +3.4 %) and their
# levels are close (median 1,942 to 2,582 kWh). KO grows a little more (+10 %
# against +5 to +7 %) and BR has more level changes (83 % against 67 to 74 %),
# but neither difference is significant with about 30 sites per suffix
# (Kruskal-Wallis p = 0.35, Fisher p = 0.37). The suffix does not call for its
# own month effect.

# ---- The two big sites -----------------------------------------------------
big_bills <- bills[bills$site %in% big, ]
p <- ggplot(big_bills, aes(month, kwh_per_day)) +
  geom_line(colour = "firebrick") + geom_point(colour = "firebrick", size = 0.8) +
  facet_wrap(~site, ncol = 1, scales = "free_y") +
  scale_y_continuous(labels = scales::label_number(big.mark = " ")) +
  scale_x_date(date_breaks = "6 months", date_labels = "%m/%Y") +
  labs(x = NULL, y = "kWh par jour")
save_fig(p, "big_sites")
# Rolling centre: 0698BR grows fast, and a site-year median would turn its
# growth into a January-to-December slope.
big_profile <- common_profile(
  seasonal_deviation(big_bills, "kwh_per_day", centre = "rolling"), by = "site")
# Conclusion: 0698BR has a strong seasonality of its own, from -11 % in
# December to +10 % in June around its rolling level (23 % between trough and
# peak), on top of +30 % growth from 2019 to 2021. 0998BR has no clear
# seasonality (within -3.4 % to +3 %) and drops about 7 % from May 2021. Neither
# follows the common profile, so they keep their own seasonality.

# ---- Naive forecast: raw kWh against kWh per day ---------------------------
# Stops before the test period of the backtest, which is read only at step 11.
naive_months <- seq(add_months(min(bills$month), 1), add_months(config$test_start, -1),
                    by = "month")
naive <- naive_errors(bills, naive_months)
naive_tab <- data.frame(
  forecast = c("naive_raw", "naive_per_day"),
  median_rel_error = c(median(naive$err_raw), median(naive$err_per_day)),
  p90_rel_error = c(quantile(naive$err_raw, 0.9), quantile(naive$err_per_day, 0.9)),
  share_per_day_better = c(NA, mean(naive$err_per_day < naive$err_raw)),
  share_equal = c(NA, mean(naive$err_per_day == naive$err_raw)),
  share_per_day_worse = c(NA, mean(naive$err_per_day > naive$err_raw))
)
# Conclusion: the per-day naive divides the median relative error by three
# (1.2 % against 3.7 %, February 2019 to May 2021; 90th percentile 6.3 % against
# 10.6 %). It is better on 70 % of the bills, equal on 14 % (two months of the
# same length) and worse on 16 %: kWh per calendar day is the right unit for
# forecasting.

# ---- Saved tables and numbers for the report -------------------------------
write.csv(sites, file.path(config$tables_dir, "site_summary.csv"), row.names = FALSE)
write.csv(profile, file.path(config$tables_dir, "month_profile.csv"), row.names = FALSE)
write.csv(suffix_tab, file.path(config$tables_dir, "suffix_summary.csv"), row.names = FALSE)
write.csv(repro, file.path(config$tables_dir, "effect_reproducibility.csv"), row.names = FALSE)
write.csv(naive_tab, file.path(config$tables_dir, "naive_raw_vs_per_day.csv"), row.names = FALSE)
write.csv(suffix_profile, file.path(config$tables_dir, "suffix_profile.csv"), row.names = FALSE)
write.csv(big_profile, file.path(config$tables_dir, "big_site_profile.csv"), row.names = FALSE)

day_prof <- profile[profile$unit == "kWh par jour", ]
raw_prof <- profile[profile$unit == "kWh du mois", ]
results <- list(
  n_sites            = nrow(sites),
  n_months           = ncol(mat),
  big_sites          = big,
  big_share_kwh      = sum(sites$share[sites$big]),
  big_median_kwh     = sites$median_kwh[sites$big],
  median_site_kwh    = median(sites$median_kwh),
  small_kwh_range    = range(sites$median_kwh[!sites$big]),
  feb_ratio_raw      = raw_prof$ratio[raw_prof$moy == 2],
  feb_ratio_day      = day_prof$ratio[day_prof$moy == 2],
  summer_ratio_day   = mean(day_prof$ratio[day_prof$moy %in% 7:8]),
  profile_range_day  = range(day_prof$ratio),
  profile_years      = full_years,
  site_effect_iqr    = stats::IQR(site_eff$effect),
  site_effect_p05_p95 = stats::quantile(exp(site_eff$effect) - 1, c(0.05, 0.95)),
  repro_median_rho   = tapply(repro$rho, repro$centre, median),
  repro_range_rho    = range(repro$rho),
  repro_dec_rho      = setNames(repro$rho[repro$moy == 12 & repro$year1 == 2019],
                                repro$centre[repro$moy == 12 & repro$year1 == 2019]),
  repro_test         = repro_test,
  repro_dec_p        = dec_p,
  suffix_tests       = suffix_tests,
  share_level_change = mean(small_sites$level_change),
  n_level_change     = sum(small_sites$level_change),
  median_growth      = median(small_sites$growth),
  naive              = naive_tab,
  naive_period       = range(naive_months),
  big_profile        = big_profile
)
saveRDS(results, file.path(config$results_dir, "descriptive.rds"))
cat("Descriptive statistics written to", config$output_dir, "\n")
