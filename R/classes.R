# Classes of behaviour of the ordinary sites (step 10). Settings in R/config.R.

#' Features of each ordinary site on the bills before `end`: median level,
#' OLS slope per year and log noise scale of y = log kWh per day minus the
#' month effect (calendar_effect), number of breaks above `min_change`.
class_features <- function(bills, big, end, seg_args, min_change) {
  b <- bills[bills$month < end & !bills$site %in% big, ]
  b <- b[order(b$site, b$month), ]
  pf <- panel_fit(b, calendar_effect(b, character()), character(), seg_args)
  large <- pf$breaks[abs(pf$breaks$change) > min_change, ]
  out <- do.call(rbind, lapply(split(pf$bills, pf$bills$site), function(d) {
    t <- seq_len(nrow(d)) / 12
    data.frame(site = d$site[1], suffix = d$suffix[1], level = stats::median(d$y),
               slope = unname(stats::coef(stats::lm(d$y ~ t))[2]),
               noise = log(d$sigma[1]), breaks = sum(large$site == d$site[1]),
               stringsAsFactors = FALSE)
  }))
  rownames(out) <- NULL
  out
}

#' Ward clustering of the standardised features (columns `cols`), number of
#' classes with the highest mean silhouette over `k_range` (smaller k on a
#' tie). Returns the class of each row, k, and the mean silhouette by k.
ward_classes <- function(features, cols, k_range) {
  x <- scale(as.matrix(features[cols]))  # (sites, features)
  d <- stats::dist(x)
  tree <- stats::hclust(d, method = "ward.D2")
  sil <- vapply(k_range, function(k) {
    mean(cluster::silhouette(stats::cutree(tree, k), d)[, "sil_width"])
  }, numeric(1))
  k <- k_range[which.max(sil)]  # which.max keeps the first, i.e. smaller, k
  list(class = unname(stats::cutree(tree, k)), k = k,
       silhouette = stats::setNames(sil, k_range), tree = tree)
}
