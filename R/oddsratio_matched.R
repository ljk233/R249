
#' @title Odd Ratio 1-1 (Matched Case-Control Study)
#' @description Return the point estimate, estimated standard error, and
#' (1-alpha)\% confidence interval for a 1-1 matched case-control study.
#'
#' @param obs Results from matched case-control study
#' @param alpha Significance level for the confidence interval.
#'
#' @export
oddsratio_matched <- function(obs, alpha = 0.05) {
  rlab <- c("1-1 matched case-control")
  clab <- c("oddsratio", "stderr", "lcb", "ucb")
  res <- matrix(nrow = 1, ncol = 4, dimnames = list(rlab, clab))
  # populate the results
  res[1, 1] <- .or_matched_est(obs)
  res[1, 2] <- .or_matched_ese(obs)
  res[1, 3] <- .or_matched_confint(obs, alpha)[1]
  res[1, 4] <- .or_matched_confint(obs, alpha)[2]
  return(res)
}


.or_matched_est <- function(obs) {
  return(obs[1, 2] / obs[2, 1])
}


.or_matched_ese <- function(obs) {
  return(sqrt(1/obs[2, 1] + 1/obs[1, 2]))
}


#' @importFrom stats qnorm
.or_matched_confint <- function(obs, alpha = 0.05) {
  or <- .or_matched_est(obs)
  z <- qnorm(1-alpha/2)
  ese <- .or_matched_ese(obs)
  c(or*exp(-z*ese), or*exp(z*ese))
}
