
#' @title Risk ratio
#'
#' @description Return the point estimate, estimated standard error, and
#' (1-alpha)\% confidence interval for the risk ratio of each exposure level.
#'
#' @param obs An r by 2 matrix containing the integer counts of the results.
#' @param alpha Significance level for the confidence interval.
#'
#' @details The reference exposure is expected in the first row, and either
#' "no disease" or cases in the first column.
#' See examples.
#'
#' @return An r by 4 matrix, one row for each exposure level
#' @export
#'
#' @examples
#' ## get results
#' res <- c(201, 76, 215, 327)
#' ## set names, optional
#' nam <- c("-", "+")
#' ## declare matrix
#' obs <- matrix(res, nrow = 2, byrow = TRUE, dimnames = list(nam, nam))
#' ## get risk ratio
#' riskratio(obs)
#' ## get risk ratio, use 1% sig level for conf ints
#' riskratio(obs, alpha = 0.01)
riskratio <- function(obs, alpha = 0.05) {
  # construct the result matrix
  nr <- nrow(obs)
  rlab <- dimnames(obs)[[1]]
  clab <- c("riskratio", "stderr", "lcb", "ucb")
  res <- matrix(nrow = nr, ncol = 4, dimnames = list(rlab, clab))
  # populate the results
  for (i in 2:nr) {
    # concat reference and exposure i
    sub <- obs[c(1, i),]
    # add results to matrix
    res[i, 1] <- .rr_est(sub)
    res[i, 2] <- .rr_ese(sub)
    res[i, 3] <- .rr_confint(sub, alpha)[1]
    res[i, 4] <- .rr_confint(sub, alpha)[2]
  }
  return(res)
}


.rr_est <- function(obs) {
  (obs[2, 2] * rowSums(obs)[1]) / (obs[1, 2] * rowSums(obs)[2])
}


.rr_ese <- function(obs) {
  sqrt((1/obs[2, 2] - 1/rowSums(obs)[2]) + (1/obs[1, 2] - 1/rowSums(obs)[1]))
}


.rr_confint <- function(obs, alpha = 0.05) {
  rr <- .rr_est(obs)
  z <- qnorm(1-alpha/2)
  ese <- .rr_ese(obs)
  c(rr*exp(-z*ese), rr*exp(z*ese))
}
