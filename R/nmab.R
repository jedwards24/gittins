#' Function arguments
#'
#' @keywords internal
#'
#' @name nmab_args
#'
#' @param Sigma Value of Sigma for the arm.
#' @param n Value of n for the arm.
#' @param gamma Numeric in (0, 1); discount factor.
#' @param tau Observation precision.
#' @param tol Absolute accuracy required.
#' @param N Time horizon used.
#' @param xi Value of xi (entent of dynamic program state space).
#' @param delta Value of delta (fineness of discretisation in the dynamic program).
NULL

#' Function arguments for value functions
#'
#' @keywords internal
#'
#' @name nmab_v_args
#'
#' @param lambda Reward from the known arm
#' @param mu mMan of reward belief for the unknown arm
#' @param n Value of n for the unknown arm
NULL

#' Calculate the Gittins index for multiple arms (normal rewards)
#'
#' Assumes mu = 0.
#'
#' @param n_range Numeric vector giving values of n (must be ascending).
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
#' The initial interval for calibration are as follows:
#' * For lower bound, use `lb` if supplied else use KGI.
#' * For upper bound, use `ub` if supplied else use GI+.
#'
#' @inheritParams nmab_args
#' @param lb Optional lower bound for GI.
#' @param ub Optional upper bound for GI.
#'
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

#' Value calculation for the one-armed bandit with Normal rewards.
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
