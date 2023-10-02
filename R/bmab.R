#' Function arguments
#'
#' @keywords internal
#' @name bmab_args
#'
#' @param Sigma Numeric > 0. Value of Sigma for the arm.
#' @param n Numeric > Sigma > 0. Value of n for the arm.
#' @param gamma Numeric in (0, 1). Reward discount factor.
#' @param N Integer >= 2. Time horizon used in the calculation..
#' @param tol Numeric > 0. Absolute accuracy required.
NULL

#' Function arguments for value functions
#'
#' @keywords internal
#' @name bmab_v_args
#'
#' @param lambda Reward from the known arm.
#' @param Sigma Value of Sigma for the unknown arm.
#' @param n Value of n for the unknown arm.
NULL

#' Calculate Gittins indices for multiple starting states (Bernoulli rewards)
#'
#' Runs `bmab_gi()` multiple times for a range of arm starting states. The starting states used are
#' determined by the starting state arguments and `num_actions` (see details). Other arguments are
#' passed to `bmab_gi()`. See `?bmab_gi` for argument and calculation details.
#'
#' @details
#' The starting states are either given in alpha/beta or Sigma/n form (strictly one or the other).
#' These are alternative ways of doing the same thing since  `alpha = Sigma` and  `beta = n - Sigma`.
#'
#' The GI are calculated for arms with the following states:
#' * `alpha = alpha_start : (alpha_start_start + num_actions - 1)`,
#' * `beta = beta_start : (beta_start + num_actions - 1)`,
#' * where `alpha + beta <= num_actions + alpha_start + beta_start - 1`.
#'
#' @inheritParams bmab_args
#' @param alpha_start Lowest value of alpha for the arms.
#' @param beta_start Lowest value of beta for the arms.
#' @param Sigma_start Lowest value of Sigma for the arms.
#' @param n_start Lowest value of n for the arms.
#' @param num_actions Determines the number of states GI are calculated for.
#'
#' @return A data frame of starting states and GI values. Extra information is stored in attributes
#' (see examples):
#' * `params`: the parameters used.
#' * `gi_matrix`: the GI values in a matrix (beta x alpha).
#' * `gi_matrix_ns`: the GI values in a matrix (n x Sigma).
#'
#' @examples
#' b1 <- bmab_gi_multiple(1, 1, num_actions = 4, gamma = 0.8, N = 20)
#' b2 <- bmab_gi_multiple(num_actions = 4, gamma = 0.8, N = 20, Sigma_start = 1, n_start = 2)
#' identical(b1, b2) # TRUE
#' # View attributes
#' attr(b1, "params")
#' attr(b1, "gi_matrix")
#' attr(b1, "gi_matrix_ns")
#'
#' @export
bmab_gi_multiple <- function(alpha_start = NULL, beta_start = NULL, num_actions, gamma, N, tol = 5e-4,
                             Sigma_start = NULL, n_start = NULL){
  if (!check_start_args(alpha_start = alpha_start, beta_start = beta_start,
                        Sigma_start = Sigma_start, n_start = n_start)){
    stop("Incorrect number of starting state arguments.\n",
         "You must supply `alpha_start` and `beta_start` OR `Sigma_start` and `n_start`.", call. = FALSE)
  }
  if (!is.null(Sigma_start)){
    if (Sigma_start >= n_start){
      stop("`n_start` must be greater than `Sigma_start`.", call. = FALSE)
    }
    alpha_start <- Sigma_start
    beta_start <- n_start - Sigma_start
  }
  check_numeric(alpha_start, "alpha_start", 0)
  check_numeric(beta_start, "beta_start", 0)
  check_integerish(num_actions, "num_actions", 2L)

  n_states <- as.integer(round(0.5 * num_actions * (num_actions + 1), 0))
  alpha_vec <- beta_vec <- integer(n_states)
  alpha_range <- alpha_start:(alpha_start + num_actions - 1)
  beta_range <- beta_start:(beta_start + num_actions - 1)
  gi_mat <- matrix(rep(NA, num_actions * num_actions), nrow = num_actions, ncol = num_actions,
                   dimnames = list(beta = beta_range, alpha = alpha_range))
  mu <- alpha_start / (alpha_start + beta_range)
  lb_vec <- bmab_kgi(alpha_start, alpha_start + beta_range, gamma)
  cat("Calculating GI values for", n_states, "states\n")
  index <- 1
  pb <- txtProgressBar(min = 0, max = num_actions, style = 3)
  for (a in 1:num_actions){
    betas <- beta_range[1:(num_actions - a + 1)]
    len_betas <- length(betas)
    end_index <- index + len_betas - 1
    alpha_vec[index:end_index] <- rep(alpha_range[a], len_betas)
    beta_vec[index:end_index] <- betas
    index <- end_index + 1
    ub <- 1
    for (b in seq_along(betas)){
      gi_mat[b, a] <- bmab_gi_ab(alpha_range[a], beta_range[b], gamma, N, tol, lb = lb_vec[b], ub = ub)
      ub <- gi_mat[b, a]
    }
    lb_vec <- gi_mat[, a]
    setTxtProgressBar(pb, a)
  }
  close(pb)
  gi_mat_ns <- gi_mat
  for(i in 2:num_actions){
    gi_mat_ns[, i] <- c(rep(NA, i - 1), gi_mat[1:(num_actions - i + 1), i])
  }
  dimnames(gi_mat_ns) <- list(n = beta_start + alpha_range, Sigma = alpha_range)
  n_vec <- alpha_vec + beta_vec
  df <- data.frame(alpha = alpha_vec, beta = beta_vec, gi = as.numeric(na.omit(as.vector(gi_mat))),
                       Sigma = alpha_vec, n = n_vec, stage = n_vec - min(n_vec))
  attr(df, "params") <- list(alpha_start = alpha_start, beta_start = beta_start,
                             Sigma_start = alpha_start, n_start = alpha_start + beta_start,
                             gamma = gamma, N = N, num_actions = num_actions, tol = tol)
  attr(df, "gi_matrix") <- gi_mat
  attr(df, "gi_matrix_ns") <- gi_mat_ns
  df
}

#' Calculate the Gittins index for a single arm (Bernoulli rewards)
#'
#' @description
#' The two versions of this function differ only in their state parameters:
#' * `bmab_gi()` uses `Sigma` and `n` parameters.
#' * `bmab_gi_ab()` uses `alpha` and `beta` parameters
#'
#' These are related by  `alpha = Sigma` (Bayesian number of successes) and  `beta = n - Sigma`
#' (Bayesian number of failures). Then `n` is the Bayesian number of observations. Together
#' with `gamma` (the discount factor for rewards), these define the problem. The remaining arguments
#' are settings for the calculation only (see details).
#'
#' @details
#' The problem has an infinite time horizon, but the dynamic program used to calculate the GI has a finite
#' horizon `N`. For sufficiently large `N`, the calculation can be arbitrarily accurate. In practice,
#' a fairly low value of `N` works well unless `gamma` is close to 1.
#'
#' The `lb` and `ub` arguments can be used to provide a starting interval for calibration if desired.
#' However, for normal use this is not needed as they will be calculated internally if not supplied.
#' So the initial interval is determined as follows:
#' * For lower bound, use `lb` if supplied else use KGI.
#' * For upper bound, use `ub` if supplied else use GI+.
#'
#' @inheritParams bmab_args
#' @param lb Optional lower bound for GI.
#' @param ub Optional upper bound for GI.
#' @seealso For a link to the accompanying paper see [gittins-package].
#'
#' @return A single Gittins index
#' @export
bmab_gi <- function(Sigma, n, gamma, N, tol = 5e-4, lb = NA, ub = NA){
  if (is.na(lb)){
    lb <- bmab_kgi(Sigma, n, gamma)
  }
  if (is.na(ub)){
    ub <- bmab_giplus(Sigma, n, gamma, tol, upper = TRUE)
  }
  check_numeric(Sigma, "Sigma", 0)
  check_numeric(n, "n", 0)
  if (Sigma >= n){
    stop("`n` must be greater than `Sigma`.", call. = FALSE)
  }
  check_numeric(gamma, "gamma", 0, 1)
  check_integerish(N, "N", 2L)
  check_numeric(tol, "tol", 0)
  check_numeric(lb, "lb")
  check_numeric(ub, "ub")
  mean(calibrate_arm(bmab_gi_value, lb, ub, tol, Sigma, n, gamma, N))
}

#' @param alpha Numeric > 0.Value of alpha for the arm.
#' @param beta Numeric > 0. Value of beta for the arm.
#'
#' @rdname bmab_gi
#' @export
bmab_gi_ab <- function(alpha, beta, gamma, N, tol = 5e-4, lb = NA, ub = NA){
  bmab_gi(Sigma = alpha, n = alpha + beta, gamma = gamma, N = N, tol = tol, lb = lb, ub = ub)
}

#' Calculate the GI+ index for a single arm (Bernoulli rewards)
#'
#' The GI+ index is an upper bound for the Gittins index.
#'
#' @inheritParams bmab_args
#' @param upper if TRUE, the upper end of the interval is returned, otherwise the midpoint.
#'
#' @return A GI+ index value.
#' @export
bmab_giplus <- function(Sigma, n, gamma, tol = 5e-4, upper = FALSE){
  interval <- calibrate_arm(bmab_giplus_value, lb = Sigma / n, ub = 1, tol, Sigma, n, gamma)
  if (upper){
    return(interval[2])
  }
  mean(interval)
}

#' Calculate the knowledge gradient index for a single arm (Bernoulli rewards)
#'
#' The KGI is an lower bound for the Gittins index.
#'
#' This is an exact closed form calculation. Arguments `Sigma` and `n` can be supplied as vectors in
#' which case a vector of index values will be returned.
#'
#' @inheritParams bmab_args
#'
#' @return An index value or a vector of values.
#' @export
bmab_kgi <- function(Sigma, n, gamma){
  mu <- Sigma / n
  H <- gamma / (1 - gamma)
  (mu + H * mu * (Sigma + 1) / (n + 1)) / (1 + H * mu)
}

#' Value of one-armed bandit using GI+ (Bernoulli rewards)
#'
#' @inheritParams bmab_v_args
#' @inheritParams bmab_args
#' @keywords internal
#' @return Difference in value between safe and unknown arms.
#' @export
bmab_giplus_value <- function(lambda, Sigma, n, gamma){
  mu <- Sigma / n
  H <- gamma / (1 - gamma)
  value_success <- H * (integrate(continue, lambda, 1, Sigma = Sigma, n = n)[[1]] +
                          lambda * pbeta(lambda, Sigma + 1, n - Sigma))
  value_fail <- H * lambda
  mu + mu * value_success + (1 - mu) * value_fail - lambda / (1 - gamma)
}

#' Helper function only used in bmab_giplus_value()
#' @noRd
continue <- function(x, Sigma, n){
  dbeta(x, Sigma + 1, n - Sigma) * x
}

#' Value of a one-armed bandit (Bernoulli rewards)
#'
#' @inheritParams bmab_v_args
#' @inheritParams bmab_args
#'
#' @keywords internal
#' @return Difference in value between safe and unknown arms.
#' @export
bmab_gi_value <- function(lambda, Sigma, n, gamma, N){
  h <- N + 1
  n_vec <- n:(n + N)
  s_vec <- Sigma:(Sigma + N)
  mu <- outer(s_vec, n_vec, "/")
  value_mat <- matrix(nrow = h, ncol = h)
  # Values of end states
  value_mat[, h] <- pmax(mu[, h], lambda) * gamma^N / (1 - gamma)
  safe_reward <- lambda * gamma^((1:N) - 1) / (1 - gamma)
  # Run DP to get values of other states
  for (i in N:1){
    j <- i + 1
    risky_reward <- mu[1:i, i] * (gamma^(i - 1) + value_mat[2:j, j]) +
      (1 - mu[1:i, i]) * value_mat[1:i, j]
    value_mat[1:i, i] <- pmax(risky_reward, safe_reward[i])
  }
  return(value_mat[1, 1] - lambda / (1 - gamma))
}
