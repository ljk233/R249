
#' @title Odds
#' @description Return the odds and log odds of an epidemiological study.
#'
#' @param obs Results from an epidemiological study.
#'
#' @export
odds <- function(obs) {
  # construct the result matrix
  nr <- nrow(obs)
  rlab <- dimnames(obs)[[1]]
  clab <- c("odds", "log(odds)")
  res <- matrix(nrow = nr, ncol = 2, dimnames = list(rlab, clab))
  # populate the results
  for (i in 2:nr) {
    res[i, 1] <- obs[i, 2] / obs[i, 1]
    res[i, 2] <- log(obs[i, 2] / obs[i, 1])
  }
  return(res)
}
