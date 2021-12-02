
#' @title Odds Ratio
#'
#' @description Return the point estimate, estimated standard error, and
#' (1-alpha)\% confidence interval for the odds ratio of each exposure level.
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
#' ## get odds ratio
#' oddsratio(obs)
#' ## get odds ratio, use 1% sig level for conf ints
#' oddsratio(obs, alpha = 0.01)
oddsratio <- function(obs, alpha = 0.05) {
  # construct the result matrix
  nr <- nrow(obs)
  rlab <- dimnames(obs)[[1]]
  clab <- c("oddsratio", "stderr", "lcb", "ucb")
  res <- matrix(nrow = nr, ncol = 4, dimnames = list(rlab, clab))
  # populate the results
  for (i in 2:nr) {
    # concat reference and exposure i
    sub <- obs[c(1, i),]
    # add results to matrix
    res[i, 1] <- .or_est(sub)
    res[i, 2] <- .or_ese(sub)
    res[i, 3] <- .or_confint(sub, alpha)[1]
    res[i, 4] <- .or_confint(sub, alpha)[2]
  }
  return(res)
}


.or_est <- function(obs) {
  (obs[1, 1] * obs[2, 2]) / (obs[1, 2] * obs[2, 1])
}


.or_ese <- function(obs) {
  sqrt(1/obs[1, 1] + 1/obs[1, 2] + 1/obs[2, 1] + 1/obs[2, 2])
}


#' @importFrom stats qnorm
.or_confint <- function(obs, alpha = 0.05) {
  or <- .or_est(obs)
  z <- qnorm(1-alpha/2)
  ese <- .or_ese(obs)
  c(or*exp(-z*ese), or*exp(z*ese))
}
