#' Performs iterative reduction of an interval using calibration with value function FUN.
#'
#' @param FUN a value function
#' @param lb initial lower end of interval
#' @param ub initial upper end of interval
#' @param tol absolute accuracy required
#'
#' @return A vector of lower and upper end of an interval containing the true value
#'
#' @export
#'
calibrate_arm <- function(FUN, lb, ub, tol, ...){
  while ((ub - lb) > tol){
    lambda <- (lb + ub) / 2
    if (FUN(lambda, ...) > 0){
      lb <- lambda
    }else{
      ub <- lambda
    }
  }
  c(lb, ub)
}
