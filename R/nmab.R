#' Function arguments
#'
#' @keywords internal
#' @name nmab_args
#'
#' @param Sigma Numeric. Value of Sigma for the arm (Bayesian reward).
#' @param n Numeric > 0. Value of n for the arm (Bayesian number of observation).
#' @param gamma Numeric in (0, 1). Reward discount factor.
#' @param tau Numeric > 0. Observation precision.
#' @param tol Numeric > 0. Absolute accuracy required.
#' @param N Integer >= 2. Time horizon used.
#' @param xi Numeric > 0. Value of xi (entent of dynamic program state space).
#' @param delta Numeric > 0. Value of delta (fineness of discretisation in the dynamic program).
NULL

#' Function arguments for value functions
#'
#' @keywords internal
#' @name nmab_v_args
#'
#' @param lambda Reward from the known arm
#' @param mu Mean of reward belief for the unknown arm
#' @param n Numeric > 0. Value of n for the unknown arm
NULL

#' Calculate Gittins indices for multiple arms (normal rewards)
#'
#' This assumes mu = 0 as the GI for other values can then be derived by simple addition. See
#' `?nmab_gi` for more detail on the problem and computation parameters.
#'
#' @param n_range Numeric vector giving values of n (all greater than 0).
#' @inheritParams nmab_args
#'
#' @return A data frame of GI vales with a row for each n in `n_range`. The parameters used are
#' attached as an attribute `params`.
#'
#' @examples
#' n1 <- nmab_gi_multiple(1:4, gamma = 0.9, tau = 1, N = 30, xi = 3, delta = 0.02, tol = 5e-4)
#' n1
#' attr(n1, "params")
#' @export
nmab_gi_multiple <- function(n_range, gamma, tau, N, xi, delta, tol = 5e-4){
  if (!is.numeric(n_range)){
    stop("`n_range` must be numeric.", call. = FALSE)
  }
  if (any(n_range <= 0)){
    stop("All values in `n_range` must be strictly greater than zero.", call. = FALSE)
  }
  check_numeric(gamma, "gamma", 0, 1)

  n_range <- sort(n_range)
  nn <- length(n_range)
  gi_vec <- numeric(nn)
  ubbl <- gamma / (1 - gamma) / sqrt(n_range[1])
  ub <- nmab_giplus(0, n_range[1], gamma, tol, ub = ubbl, upper = TRUE)
  cat("Calculating GI values for", nn, "states (may be slow)\n")
  pb <- txtProgressBar(min = 0, max = nn, style = 3)
  for (i in 1:nn){
    gi_vec[i] <- nmab_gi(0, n_range[i], gamma, tau, N, xi, delta, tol = tol, ub = ub)
    ub <- gi_vec[i]
    setTxtProgressBar(pb, i)
  }
  close(pb)
  df <- data.frame(n = n_range, gi = gi_vec)
  attr(df, "params") <- list(n_range = n_range, gamma = gamma, tau = tau,
                             N = N, xi = xi, delta = delta, tol = tol)
  df
}

#' Calculate the Gittins index for a single arm (normal rewards)
#'
#' The problem state is given by the `Sigma`, `n`, `gamma`, and `tau` arguments. The
#' remaining arguments affect how the calculation is done and are chosen based on accuracy
#' and speed requirements (see details).
#'
#' @details
#' The problem has an infinite continuous state space, but the calculation can only be done for a
#' finite discrete area. `N` and `xi` control the finite extent of the state space used (one in
#' each direction). These have a diminishing effect on accuracy and just need to be sufficiently
#' large, beyond which there will be minimal improvement. Each parameter has a linear effect on the state
#' space size and therefore algorithm speed. Larger `gamma` or smaller `tau` requires
#' a larger `N`, but depends only on the least favourable of the two. Often `N` can be smaller
#' than for similar BMAB problems. The required `xi` is more robust to changes in the problem setting
#' as adjusts to `tau` somewhat. Therefore, standard values can be used across problems with less
#' disadvantage. `xi` is based on standard deviations so 3 or 4 will often be large enough. Using `N` or
#' `xi` values that are too low gives an underestimate of the Gittins index.
#'
#' The discretisation parameter `delta` is the main limiter to accuracy as it has a non-linear
#' effect on computation time. Accuracy increases as `delta` gets smaller, with inaccuracies leading
#' to an overestimate of the Gittins index.
#'
#' The `lb` and `ub` arguments can be used to provide a starting interval for calibration if desired.
#' However, for normal use this is not needed as they will be calculated internally if not supplied.
#' So the initial interval is determined as follows:
#' * For lower bound, use `lb` if supplied else use KGI.
#' * For upper bound, use `ub` if supplied else use GI+.
#'
#' @inheritParams nmab_args
#' @param lb Optional lower bound for GI.
#' @param ub Optional upper bound for GI.
#' @seealso For a link to the accompanying paper see [gittins-package].
#' @return A vector of GI values.
#'
#' @export
nmab_gi <- function(Sigma, n, gamma, tau, N, xi, delta, tol = 5e-4, lb = NA, ub = NA){
  if (is.na(lb)){
    lb <- nmab_kgi(0, n, gamma, tau, tol, ub, lower = TRUE)
  }
  if (is.na(ub)){
    ub <- nmab_giplus(0, n, gamma, tol, ub, upper = TRUE)
  }
  check_numeric(Sigma, "Sigma")
  check_numeric(n, "n", 0)
  check_numeric(gamma, "gamma", 0, 1)
  check_numeric(tau, "tau", 0)
  check_integerish(N, "N", 2L)
  check_numeric(xi, "xi", 0)
  check_numeric(delta, "delta", 0)
  check_numeric(tol, "tol", 0)
  check_numeric(lb, "lb")
  check_numeric(ub, "ub")

  interval <- Sigma / n + calibrate_arm(nmab_gi_value, lb, ub, tol, n, gamma, tau, N, xi, delta)
  mean(interval)
}

#' Calculate the GI+ index for a single arm (normal rewards)
#'
#' The index is an upper bound for GI.
#'
#' @param upper If TRUE, the upper end of the interval is returned, otherwise the midpoint.
#' @param ub Optional upper bound for the index.
#' @inheritParams nmab_args
#'
#' @return A vector of GI values.
#' @export
nmab_giplus <- function(Sigma, n, gamma, tol = 5e-4, ub = NA, upper = FALSE){
  if(is.na(ub)){ub <- gamma / (1 - gamma) / sqrt(n)}
  interval <- Sigma / n + calibrate_arm(nmab_giplus_value, lb = 0, ub = ub, tol, mu = 0, n, gamma)
  if (upper){
    return(interval[2])
  }
  mean(interval)
}

#' Calculate the knowledge gradient index for a single arm (normal rewards)
#'
#' The index is an lower bound for GI.
#'
#' @inheritParams nmab_args
#' @param lower If TRUE, the lower end of the interval is returned, otherwise the midpoint.
#' @param ub Optional upper bound for the index.
#'
#' @return A vector of GI values.
#' @export
nmab_kgi <- function(Sigma, n, gamma, tau, tol = 5e-4, ub = NA, lower = FALSE){
  if(is.na(ub)){ub <- gamma / (1 - gamma) / sqrt(n)}
  interval <- Sigma / n + calibrate_arm(nmab_kgi_value, lb = 0, ub, tol, mu = 0, n, gamma, tau)
  if (lower){
    return(interval[1])
  }
  mean(interval)
}

#' Value of one-armed bandit using GI+ (normal rewards)
#'
#' @inheritParams nmab_v_args
#' @inheritParams nmab_args
#' @keywords internal
#'
#' @return Difference in value between safe and unknown arms.
#' @export
nmab_giplus_value <- function(lambda, mu, n, gamma){
  sd <- sqrt(1 / n)
  mu + gamma * dnorm(lambda / sd) * sd / (1 - gamma * pnorm(lambda / sd)) - lambda
}

#' Value of one-armed bandit using KGI (normal rewards)
#'
#' Calculation is as given in Ryzhov, Powell & Frazier (2012), simplified since `lambda >= mu`
#' always when called here.
#'
#' @inheritParams nmab_v_args
#' @inheritParams nmab_args
#' @keywords internal
#'
#' @return Difference in value between safe and unknown arms.
#' @export
nmab_kgi_value <- function(lambda, mu, n, gamma, tau){
  sigt <- sqrt(1 / n - 1 / (n + tau))
  z <- (mu - lambda) / sigt
  v <- sigt * (z * pnorm(z) + dnorm(z))
  v * gamma / (1 - gamma) + mu - lambda
}

#' Reward of the risky arm in a one-armed bandit process
#'
#' Helper function only used in nmab_gi_value.
#' @noRd
nmab_risky_reward <- function(mu, y_lo_scaled, y_hi_scaled, tn_scaled, tau, s, value_vec, discount){
  yhi <- y_hi_scaled - mu * tn_scaled
  ylo <- y_lo_scaled - mu * tn_scaled
  p <- pnorm(yhi, mu, s) - pnorm(ylo, mu, s)
  discount * mu + (pnorm(ylo[1], mu, s) * value_vec[1] +
                     (1 - pnorm(yhi[length(yhi)], mu, s)) * value_vec[length(value_vec)] +
                     sum(p * value_vec))
}

#' Value calculation for the one-armed bandit with Normal rewards
#'
#' Assumes `Sigma = mu = 0`.
#'
#' The `extra_xi` argument was a later addition to the algorithm, not included in the paper, which
#' improves accuracy at low computational cost.
#'
#' Normally, states outside the width of the state space are ignored (taken to have a value of zero).
#' This saves computation for states that are unlikely to be visited. However, the calculation can
#' be improved with relatively little computation by giving some of these states a value using
#' their mean reward only (no further learning). Although this is an approximation it will always
#' be more accurate than using zero. So there are two blacks of states: the original states are
#' within xi standard deviations and are calculated in detail using dynamic programming; and the
#' new states within xi + extra_xi standard deviations. I have found `extra_xi = 1` works well and
#' have set this as the default. This is value that should be used unless doing research on its effect.
#'
#' @inheritParams nmab_v_args
#' @inheritParams nmab_args
#' @param extra_xi Extend xi using a fast approximation. See details
#' @keywords internal
#'
#' @return Difference in value between safe and unknown arms.
#' @export
nmab_gi_value <- function(lambda, n, gamma, tau, N, xi, delta, extra_xi = 1){
  h <- N + 1
  delta <- delta / sqrt(n) # adjust delta so number of states is constant with n
  mu_range <- seq(0, (xi + extra_xi) * sqrt(1 / n), by = delta)
  lr <- length(mu_range)
  lr2 <- length(seq(0, (xi) * sqrt(1 / n), by = delta))
  value <- matrix(nrow = h, ncol = lr)
  # Value of end states (at stage N)
  value[h, ] <- pmax(mu_range, lambda) * gamma^N / (1 - gamma)
  rr <- gamma^(0:N) / (1 - gamma)
  value[, (lr2:lr)] <- matrix(rr, ncol = 1) %*% matrix(mu_range[lr2:lr], nrow = 1)
  lo <- mu_range - delta / 2 # lower end of discrete approximation to mu
  hi <- mu_range + delta / 2
  for (j in N:2){
    t <- j - 1
    tn <- n + tau * (j - 1)
    #the next 3 variables are used for speed-up only
    y_hi_scaled <- hi * (tn + tau) / tau
    y_lo_scaled <- lo * (tn + tau) / tau
    tn_scaled <- tn / tau
    s <- sqrt(1 / tn + 1 / tau) #sd of y
    discount <- gamma^t
    safe_reward <- lambda * discount / ( 1 - gamma)
    value_vec <- value[j + 1, ]
    for (i in (lr2 - 1):1){
      risky_reward <- nmab_risky_reward(mu_range[i], y_lo_scaled, y_hi_scaled, tn_scaled,
                                        tau, s, value_vec, discount)
      if (risky_reward > safe_reward){
        value[j, i] <- risky_reward
      }else{
        value[j, 1:i] <- safe_reward
        break
      }
    }
    # If risky arm preferred in [j, 1] then it will be preferred in starting state
    if (value[j, 1] > safe_reward){
      return(value[j, 1] - safe_reward)
    }
  }
  # Value of risky arm in starting state at time 0
  s <- sqrt(1 / n + 1 / tau)
  value_vec = value[2, ]
  risky_reward <- nmab_risky_reward(mu_range[1],
                                    y_lo_scaled = lo * (n + tau) / tau,
                                    y_hi_scaled = hi * (n + tau) / tau,
                                    tn_scaled = n / tau,
                                    tau, s, value_vec, discount = 1)
  return(risky_reward - lambda / (1 - gamma))
}
