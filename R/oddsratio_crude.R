
#' @title Crude Odds Ratio
#'
#' @description Return the pooled unweighted odds ratio for a stratified analysis.
#'
#' @param obs An r by 2 by z matrix.
#' @param alpha Significance level for the confidence interval.
#'
#' @return An r by 4 matrix, one row for each exposure level.
#' @export
oddsratio_crude <- function(obs, alpha = 0.05){
  nr <- nrow(obs)
  nz <- dim(obs)[3]
  zeroes <- rep(c(0), nr*2)
  res <- matrix(zeroes, nrow = nr, ncol = 2)
  # aggregate obs over levels
  for (z in 1:nz) {
    for (r in 1:nr) {
      for (c in 1:2) {
        res[r, c] <- res[r, c] + obs[r, c, z]
      }
    }
  }
  return(oddsratio(res, alpha))
}
