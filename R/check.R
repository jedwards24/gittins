#' Check if a vector is a scalar numeric strictly within given bounds
#'
#' Returns `NULL` if okay. Errors otherwise.
#' @param x Object to test.
#' @param name String to use as argument name in error message. Will be enclosed in back ticks.
#' @param min Lower bound for `x`. `x > min`.
#' @param max Upper bound for `x`. `x < max`.
#' @noRd
check_numeric <- function(x, name, min = -Inf, max = Inf, int = FALSE) {
  name <- paste0("`", name, "`")
  if (length(x) > 1 || !is.numeric(x)){
    stop(name, " must be a length 1 numeric.", call. = FALSE)
  }
  if (x <= min){
    stop(name, " must be greater than ", min, ".", call. = FALSE)
  }
  if (x >= max){
    stop(name, " must be less than ", max, ".", call. = FALSE)
  }
  return()
}

#' Check if a vector is a scalar integer like numeric within given bounds
#'
#' Unlike `check_numeric()` this allows equality with bounds.
#' Returns `NULL` if okay. Errors otherwise.
#' @param x Object to test.
#' @param name String to use as argument name in error message. Will be enclosed in back ticks.
#' @param min Lower bound for `x`. `x >= min`.
#' @param max Upper bound for `x`. `x <= max`.
#' @param int Logical. Should `x` be integer like?
#' @noRd
check_integerish <- function(x, name, min = -Inf, max = Inf) {
  name <- paste0("`", name, "`")
  if (length(x) > 1 || !is.numeric(x)){
    stop(name, " must be a length 1 numeric.", call. = FALSE)
  }
  if (x < min){
    stop(name, " must be at least ", min, ".", call. = FALSE)
  }
  if (x > max){
    stop(name, " must be no greater than ", max, ".", call. = FALSE)
  }
  if (!round(x, 0) == x){
    stop(name, " must be a whole number.", call. = FALSE)
  }
  return()
}

#' Check if the correct starting state arguments have values
#'
#' Helper for `bmab_gi_multiple()`. Returns `TRUE` if arguments are okay, `FALSE` otherwise.
#' @noRd
check_start_args <- function(alpha_start = NULL, beta_start = NULL, Sigma_start = NULL, n_start = NULL) {
  if (is.null(alpha_start) + is.null(beta_start) == 1) return(FALSE)
  if (is.null(Sigma_start) + is.null(n_start) == 1) return(FALSE)
  return(xor(!is.null(alpha_start) && !is.null(beta_start),
             !is.null(Sigma_start) && !is.null(n_start)))
}
