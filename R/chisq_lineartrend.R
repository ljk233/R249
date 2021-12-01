
#' @title Chi-Squared Test for No Linear Trend (Dose-Response Analysis)
#' @description Return the test statistic and p-value for a chi-squared test for
#' no linear trend.
#' @references https://online.stat.psu.edu/stat504/book/export/html/710
#' @param tbl Results from a dose-response study.
#' @importFrom stats pchisq
#' @export
chisq_lineartrend <- function(tbl) {
  # gather results
  rscore <- .midranks(tbl)
  cscore <- c(1, 2)
  rbar <- sum(margin.table(tbl,1) * rscore)/sum(tbl)
  cbar <- sum(margin.table(tbl,2) * cscore)/sum(tbl)
  rdif <- rscore - rbar
  cdif <- cscore - cbar
  sdr <- sum(margin.table(tbl,1) * (rdif^2))
  sdc <- sum(margin.table(tbl,2) * (cdif^2))
  sdrc <- sum(t(tbl * rdif) * cdif)
  r <- sdrc/(sqrt(sdr * sdc))
  # calculate M2, p-value
  M2 <- (sum(tbl) - 1) * r^2
  pval <- pchisq(M2, 1, lower.tail = FALSE)
  return(
    matrix(c(M2, pval), ncol = 2, dimnames = list(c("result"), c("chisq", "pval")))
  )
}


.midranks <- function(obs) {
  ranks <- rep(0, times = dim(obs)[1])
  rowsums <- rep(0, times = dim(obs)[1])
  for (i in 1:dim(obs)[1]) {
    rowsums[i] <- obs[i, 1] + obs[i, 2]
  }
  counter <- 0
  for (i in 1:length(rowsums)) {
    ranks[i] <- ranks[i] + (counter + ((1 + rowsums[i]) / 2))
    counter <- counter + rowsums[i]
  }
  return(ranks)
}
