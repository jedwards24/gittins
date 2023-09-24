#' Calibrate an unknown arm against a known arm
#'
#' Performs iterative reduction of an interval using calibration with value function FUN.
#'
#' @param FUN A value function.
#' @param lb Initial lower end of interval.
#' @param ub Initial upper end of interval.
#' @param tol Absolute accuracy required.
#' @param ... Other arguments passed to `FUN`.
#'
#' @return A vector of lower and upper end of an interval containing the true value.
#' @export
calibrate_arm <- function(FUN, lb, ub, tol, ...){
  while ((ub - lb) > tol){
    lambda <- lb + (ub - lb) / 2
    if (FUN(lambda, ...) > 0){
      lb <- lambda
    }else{
      ub <- lambda
    }
  }
  c(lb, ub)
}
