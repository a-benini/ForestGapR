#' Forest Canopy Gap-size Frequency Distribution
#'
#' @description This function quantifies the forest canopy gap size-frequency distribution and estimates the power-law exponent (\ifelse{latex}{\out{$\lambda$}}{\ifelse{html}{\out{&lambda;}}{lambda}}) from the Zeta distribution.
#'
#' @usage GapSizeFDist(gaps_stats,method,...)
#'
#' @param gaps_stats A \code{data.frame} containing basic statistics of forest gaps. Output of [GapStats()] function.
#' @param method method computing the \ifelse{latex}{\out{$\lambda$}}{\ifelse{html}{\out{&lambda;}}{lambda}}: \code{Hanel_2017} (default) described in Hanel et al. (2017) or \code{Asner_2013} described in Asner et al. (2013).
#' @param ... Supplementary parameters for [graphics::plot()].
#' @return A log-log plot of the gap-size frequency distribution and a list containing: i) \ifelse{latex}{\out{$\lambda$}}{\ifelse{html}{\out{&lambda;}}{lambda}}, ii) the gap-size frequency distribution and iii) the method used.
#' The \ifelse{latex}{\out{$\lambda$}}{\ifelse{html}{\out{&lambda;}}{lambda}} parameter is derived for the Zeta distribution using a maximum likelihood estimator. See details section.
#'
#' @references
#'
#' Hanel,R., Corominas-Murtra, B., Liu, B., Thurner, S. (2013). Fitting power-laws in empirical data with estimators that work for all exponents,
#' PloS one, vol. 12, no. 2, p. e0170920.https://doi.org/10.1371/journal.pone.0170920
#'
#' Asner, G.P., Kellner, J.R., Kennedy-Bowdoin, T., Knapp, D.E., Anderson, C. & Martin, R.E. (2013).
#' Forest canopy gap distributions in the Southern Peruvian Amazon. PLoS One, 8, e60875.https://doi.org/10.1371/journal.pone.0060875
#'
#' White, E.P, Enquist, B.J, Green, J.L. (2008) On estimating the exponent of power law frequency distributions. Ecology 89,905-912.
#' https://doi.org/10.1890/07-1288.1
#'
#' @examples
#' # Loading terra library
#' library(terra)
#'
#' # ALS-derived CHM over Adolpho Ducke Forest Reserve - Brazilian tropical forest
#' ALS_CHM_DUC <- rast(system.file("tif/ALS_CHM_DUC.tif", package = "ForestGapR"))
#'
#' # set height thresholds (e.g. 10 meters)
#' threshold <- 10
#' size <- c(1, 10^4) # m2
#'
#' # Detecting forest gaps
#' gaps_duc <- getForestGaps(chm_layer = ALS_CHM_DUC, threshold = threshold, size = size)
#'
#' # Computing basic statistics of forest gap
#' gaps_stats <- GapStats(gap_layer = gaps_duc, chm_layer = ALS_CHM_DUC)
#'
#' # Gap-size Frequency Distributions
#' GapSizeFDist(
#'   gaps_stats = gaps_stats, method = "Hanel_2017", col = "forestgreen", pch = 16, cex = 1,
#'   axes = FALSE, ylab = "Gap Frequency", xlab = as.expression(bquote("Gap Size" ~ (m^2)))
#' )
#' axis(1)
#' axis(2)
#' grid(4, 4)
#' @export
GapSizeFDist <- function(gaps_stats, method = "Hanel_2017", ...) {
  method <- match.arg(method, c("Hanel_2017", "Asner_2013"))
  if (method == "Asner_2013") {
    fit <- stats::optimize(function(data, lambda) {
      2 * sum(-log(data^-lambda / VGAM::zeta(x = lambda)))
    }, data = gaps_stats$gap_area, lower = 1.0001, upper = 20, maximum = F)
    .lambda <- fit$minimum
  }else{ # method = "Hanel_2017" / default
    a <- poweRlaw::conpl$new(gaps_stats$gap_area)
    .lambda <- as.numeric(poweRlaw::estimate_pars(a)[1]) # saving as a separate object
  }
  
  gap.size <- seq(0, max(gaps_stats$gap_area), 1)
  gap.freq <- as.numeric(table(base::cut(gaps_stats$gap_area, breaks = gap.size)))
  gap.freq_mat <- cbind(gap.size = gap.size[-1], gap.freq = gap.freq)
  
  graphics::plot(gap.freq~gap.size, data = gap.freq_mat[gap.freq_mat[, "gap.freq"] > 0, ],  log = "xy", ...)
  eqn <- bquote(lambda == .(round(.lambda, 3)) * "," ~ ~ n == .(nrow(gaps_stats)))
  graphics::legend("topright", legend = eqn, bty = "n")

  return(list(lambda = .lambda, gap.freq = gap.freq_mat, method = method))
}
