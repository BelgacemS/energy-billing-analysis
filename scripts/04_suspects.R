# Step 6: typology of the flagged bills and gaps in kWh.
# Needs the outputs of steps 4 and 5 (baseline.rds, panel.rds).
# Run from the project root: Rscript scripts/04_suspects.R
source("R/config.R")
source("R/classify.R")
suppressPackageStartupMessages(library(ggplot2))

for (d in c(config$figures_dir, config$tables_dir, config$results_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}
save_fig <- function(plot, name, width = 7, height = 4.5) {
  ggsave(file.path(config$figures_dir, paste0(name, ".png")), plot,
         width = width, height = height, dpi = 300, bg = "white")
}
theme_set(theme_minimal(base_size = 10))

panel <- readRDS(file.path(config$results_dir, "panel.rds"))
periods <- readRDS(file.path(config$results_dir, "baseline.rds"))$periods
pb <- panel$panel
big <- unique(pb$site[pb$group != "common"])  # big sites keep their own month effect

# ---- Typology ----------------------------------------------------------------
flagged <- rank_flagged(classify_bills(pb, panel$breaks, periods, big, config))
# Round trips (rule added after review, outside the step 7 calibration).
trips <- find_round_trips(panel$breaks[!panel$breaks$site %in% big, ],
                          config$plateau_min_months, config$level_shift_min)
flagged <- apply_round_trips(flagged, pb, trips, config)
# Lists of the report: probable billing errors, changes of the site to
# confirm (a shared break stays here: the per-month test is not corrected for
# the number of months tested), indeterminate; ramps are shown, not listed.
list_of <- c(irregular_period = "error", point = "error", implausible = "error",
             temporary = "error", compensated_pair = "error", round_trip = "error",
             lasting_change = "change", network_change = "change",
             indeterminate = "indeterminate", transition = "none")
flagged$list <- unname(list_of[flagged$type])
stopifnot(!anyNA(flagged$list))
small <- flagged[!flagged$big, ]
big_tab <- flagged[flagged$big, ]

type_labels <- c(irregular_period = "Période irrégulière", point = "Écart ponctuel",
                 implausible = "Valeur implausible", temporary = "Écart temporaire",
                 compensated_pair = "Paire compensée", round_trip = "Aller-retour de niveau",
                 lasting_change = "Changement durable",
                 indeterminate = "Indéterminé", transition = "Transition",
                 network_change = "Probable changement réseau")
conf_labels <- c(high = "élevée", medium = "moyenne", low = "faible")

by_type_conf <- table(type = factor(small$type, levels = names(type_labels)),
                      confidence = small$confidence)
episodes <- unique(small[!is.na(small$episode), c("episode", "type")])
# Gaps of the suspect bills (break rows apart: their gap is against the old
# level and their impact is the cumulative gap of the lasting changes).
bill_rows <- small$suspect & (!is.na(small$episode) | small$type == "irregular_period")
gap_by_direction <- tapply(small$gap_kwh[bill_rows], small$direction[bill_rows], sum)
changes <- small[small$type == "lasting_change", ]
cumul_by_direction <- tapply(changes$cumul_kwh, changes$direction, sum)

# Validated rule check: no temporary run longer than temporary_max_months.
stopifnot(longest_temporary(small) <= config$temporary_max_months)

# Shared breaks against chance (same call as inside classify_bills).
testable <- testable_months(pb$month, config$plateau_min_months)
shared_months <- shared_break_months(panel$breaks, big, testable, config$level_shift_min,
                                     config$shared_break_min_sites, config$shared_break_alpha)

# Compensated pairs, one row per episode: months, order and net gap at the
# level that validated the pair (the gaps month by month use each bill's own
# plateau and do not sum to the net).
pr <- small[small$type == "compensated_pair", ]
pairs <- do.call(rbind, lapply(split(pr, pr$episode), function(e) {
  data.frame(episode = e$episode[1], site = e$site[1], first = min(e$month), last = max(e$month),
             n_bills = nrow(e), order = e$pair_order[1], net_kwh = e$pair_net_kwh[1],
             sum_monthly_gaps = sum(e$gap_kwh))
}))


# 0089KO 2019: June to August against the level that follows (subject: the
# three bills sum to about three months of the new level).
p89 <- pb[pb$site == "0089KO" & pb$month >= as.Date("2019-06-01") &
            pb$month <= as.Date("2019-08-01"), ]
level_after <- pb$level[pb$site == "0089KO" & pb$month == as.Date("2019-09-01")]
# expected_after: the new level from June (the level that validates the pair);
# expected_own: each bill at its own plateau (old level until July).
pair_0089 <- c(billed = sum(p89$kwh),
               expected_after = sum(exp(level_after + p89$effect) * p89$month_days),
               expected_own = sum(p89$expected_kwh))
pair_0089["ratio"] <- pair_0089[["billed"]] / pair_0089[["expected_after"]]

case_sites <- c("0089KO", "1067KO", "1088KO", "1978KO", "1986KO", "0254BB", "1240BR", "1911BB")
cases <- flagged[flagged$site %in% case_sites, ]
cases <- cases[order(cases$site, cases$month), ]

# Values behind the qualitative sentences of the report.
suspect_rows <- small[small$suspect, ]
is_break <- is.na(small$episode) & small$type != "irregular_period"
returning <- small$type %in% c("point", "temporary", "implausible")
p0089_may <- pb[pb$site == "0089KO" & pb$month == as.Date("2019-05-01"), ]
facts <- list(
  n_bills_suspect   = nrow(unique(suspect_rows[c("site", "month")])),
  types_by_count    = sort(table(small$type), decreasing = TRUE),
  suspect_types_by_count = sort(table(suspect_rows$type), decreasing = TRUE),
  # months of each type for the cases of the subject
  case_months       = stats::aggregate(month ~ site + type, cases[!cases$big, ],
                                       function(m) paste(format(sort(m), "%Y-%m"), collapse = ", ")),
  breaks_1240BR     = cases[cases$site == "1240BR" & is.na(cases$episode),
                            c("month", "type", "change", "cumul_kwh", "cumul_from")],
  may2019_0089KO    = c(z = p0089_may$z, flagged = p0089_may$excursion),
  # billed / expected: bills back to their level vs break rows
  ratio_returning   = range(small$kwh[returning] / small$expected_kwh[returning]),
  ratio_breaks      = range(small$kwh[is_break] / small$expected_kwh[is_break]),
  big_edge          = big_tab[c("site", "month", "edge", "confidence")]
)

# The three lists of the report, ordinary sites. The gaps of the change list
# are taken at the old level on the first month only: counts, not kWh.
lists <- do.call(rbind, lapply(c("error", "change", "indeterminate"), function(l) {
  x <- small[small$list == l, ]
  data.frame(list = l, rows = nrow(x), bills = nrow(unique(x[c("site", "month")])),
             sites = length(unique(x$site)),
             over_kwh = sum(x$gap_kwh[x$gap_kwh > 0]), under_kwh = sum(x$gap_kwh[x$gap_kwh < 0]))
}))
listed <- small[small$list != "none", ]
list_facts <- list(
  lists = lists,
  types = table(factor(listed$type, levels = names(type_labels))),
  n_bills = nrow(unique(listed[c("site", "month")])),
  n_sites = length(unique(listed$site)),
  n_trips = nrow(trips),
  trip_sites = trips$site,
  kwh_ordinary = sum(pb$kwh[!pb$site %in% big]),
  n_bills_ordinary = sum(!pb$site %in% big)
)

# ---- Figure ------------------------------------------------------------------
counts <- as.data.frame(by_type_conf, responseName = "n")
counts$type <- factor(type_labels[as.character(counts$type)], levels = rev(type_labels))
counts$confidence <- factor(conf_labels[as.character(counts$confidence)], levels = conf_labels)
p <- ggplot(counts, aes(n, type, fill = confidence)) +
  geom_col() +
  scale_fill_manual(values = c("#08519c", "#6baed6", "#c6dbef"), name = "Confiance") +
  labs(x = "Nombre de lignes (factures signalées et ruptures)", y = NULL,
       caption = "Sites ordinaires. Transition et probable changement réseau : affichés, non suspects.") +
  theme(legend.position = "bottom")
save_fig(p, "suspects_by_type")
# Conclusion (K = 6 from step 7): on the 90 ordinary sites, 3 round trips
# (0357KO 2020, 0690BR and 2106BB 2021; 0254BB and 1037BB alternate up and
# down and stay ambiguous). Lists: 44 probable errors on 26 sites (19 point
# gaps, 10 temporary bills, 9 round-trip bills, 4 bills in 2 compensated
# pairs, 2 irregular periods, no implausible value), +9,075 and -6,307 kWh
# against the expected value; 85 changes to confirm on 57 sites (56 lasting,
# 29 in the 5 month-directions that beat chance, a test not corrected for the
# 60 month-directions tested); 14 indeterminate rows (13 bills) on 8 sites.
# 27 ramp bills are shown, not listed. The big sites have no listed row.

# ---- Saved tables and numbers for the report ---------------------------------
write.csv(small, file.path(config$tables_dir, "suspects.csv"), row.names = FALSE)
write.csv(big_tab, file.path(config$tables_dir, "suspects_big.csv"), row.names = FALSE)

results <- list(
  config = config[c("segment_K", "level_shift_min", "pair_tolerance", "implausible_ratio",
                    "shared_break_min_sites", "shared_break_alpha", "edge_months", "merge_gap_months",
                    "temporary_max_months", "high_confidence_z", "range_sigmas")],
  type_labels       = type_labels,
  conf_labels       = conf_labels,
  flagged           = flagged,
  n_rows            = nrow(small),
  n_suspect         = sum(small$suspect),
  n_sites           = length(unique(small$site)),
  n_suspect_sites   = length(unique(small$site[small$suspect])),
  by_type_conf      = by_type_conf,
  by_type           = rowSums(by_type_conf),
  by_conf_suspect   = table(small$confidence[small$suspect]),
  episodes_by_type  = table(factor(episodes$type, levels = names(type_labels))),
  gap_by_direction  = gap_by_direction,
  cumul_by_direction = cumul_by_direction,
  big               = big_tab,
  pair_0089         = pair_0089,
  cases             = cases,
  n_bills_suspect   = facts$n_bills_suspect,
  n_testable_months = length(testable),
  shared_lambda     = unique(shared_months[c("direction", "lambda")]),
  shared_months     = shared_months,  # month, direction, n, lambda, p_value, shared
  pairs             = pairs,
  facts             = facts,
  trips             = trips,
  list_facts        = list_facts
)
saveRDS(results, file.path(config$results_dir, "suspects.rds"))
cat("Suspects written to", config$output_dir, "\n")
