#' Function arguments
#'
#' @keywords internal
#'
#' @name nmab_args
#'
#' @param Sigma value of Sigma for the arm
#' @param n value of n for the arm
#' @param gamma numeric in (0, 1); discount factor
#' @param tau observation precision
#' @param tol absolute accuracy required
#' @param N time horizon used
#' @param xi value of xi (entent of dynamic programmme state space)
#' @param delta value of delta (fineness of discretisation in the dynamic programmme)
#'
NULL

#' Function arguments for value functions
#'
#' @keywords internal
#'
#' @name nmab_v_args
#'
#' @param lambda reward from the known arm
#' @param mu mean of reward belief for the unknown arm
#' @param n value of n for the unknown arm
#'
NULL

#' Calculate the Gittins index for multiple arms (normal rewards)
#'
#' Assumes mu=0.
#'
#' @param n_range numeric vector giving values of n (must be ascending)
#' @inheritParams nmab_args
#'
#' @return A vector of GI values
#'
#' @examples
#' nmab_gi_multiple(1 : 20, gamma = 0.9, tau = 1, tol = 5e-5, N = 30, xi = 3, delta = 0.02)
#'
#' @export
#'
nmab_gi_multiple <- function(n_range, gamma, tau, tol, N, xi, delta){
  nn <- length(n_range)
  gi_vec <- numeric(nn)
  ubbl <- gamma / (1 - gamma) / sqrt(n_range[1])
  ub <- nmab_giplus(0, n_range[1], gamma, tol, ub=ubbl, upper=T)
  cat("Calculating GI values for", nn, "states (may be slow)\n")
  pb <- txtProgressBar(min = 0, max = nn, style = 3)
  for (i in 1 : nn){
    gi_vec[i] <- nmab_gi(0, n_range[i], gamma, tau, tol, N, xi, delta, ub=ub)
    ub <- gi_vec[i]
    setTxtProgressBar(pb, i)
  }
  close(pb)
  gi_vec
}

#' Calculate the Gittins index for a single arm (normal rewards)
#'
#' The initial interval for calibration are as follows:
#' For lower bound, use lb if supplied else use KGI.
#' For upper bound, use ub if supplied else use GI+.
#'
#' @inheritParams nmab_args
#' @param lb=NA optional lower bound for GI
#' @param ub=NA optional upper bound for GI
#'
#' @return A vector of GI values
#'
#' @export
#'
nmab_gi <- function(Sigma, n, gamma, tau, tol, N, xi, delta, lb=NA, ub=NA){
  if (is.na(lb)){
    lb <- nmab_kgi(0, n, gamma, tau, tol, ub, lower=T)
  }
  if (is.na(ub)){
    ub <- nmab_giplus(0, n, gamma, tol, ub, upper=T)
  }
  interval <- Sigma / n + calibrate_arm(nmab_gi_value, lb, ub, tol, n, gamma, tau, N, xi, delta)
  mean(interval)
}

#' Calculate the GI+ index for a single arm (normal rewards)
#'
#' The index is an upper bound for GI.
#'
#' @param upper=F if TRUE, the upper end of the interval is returned, otherwise the midpoint
#' @inheritParams nmab_args
#'
#' @return A vector of GI values
#'
#' @export
#'
nmab_giplus <- function(Sigma, n, gamma, tol, ub = NA, upper = F){
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
#' @param lower=F if TRUE, the lower end of the interval is returned, otherwise the midpoint
#'
#' @return A vector of GI values
#'
#' @export
#'
nmab_kgi <- function(Sigma, n, gamma, tau, tol, ub = NA, lower=F){
  if(is.na(ub)){ub <- gamma / (1 - gamma) / sqrt(n)}
  interval <- Sigma / n + calibrate_arm(nmab_kgi_value, lb=0, ub, tol, mu=0, n, gamma, tau)
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
#' @return Difference in value between safe and unknown arms
#'
#' @export
#'
nmab_giplus_value <- function(lambda, mu, n, gamma){
  sd <- sqrt(1 / n)
  mu + gamma * dnorm(lambda / sd) * sd / (1 - gamma * pnorm(lambda / sd)) - lambda
}

#' Value of one-armed bandit using KGI (normal rewards)
#'
#' Calculation is as given in Ryzhov, Powell & Frazier (2012), simplified since `lambda>=mu` always when called here.
#'
#' @inheritParams nmab_v_args
#' @inheritParams nmab_args
#'
#' @return Difference in value between safe and unknown arms
#'
#' @export
#'
nmab_kgi_value <- function(lambda, mu, n, gamma, tau){
  sigt <- sqrt(1 / n - 1 / (n + tau))
  z <- (mu - lambda) / sigt
  v <- sigt * (z * pnorm(z) + dnorm(z))
  v * gamma / (1 - gamma) + mu - lambda
}


#' Reward of the risky arm in a one-armed bandit process
#'
#' Helper function only used in nmab_gi_value.
#'
#' @keywords internal
#'
nmab_risky_reward <- function(mu, y_lo_scaled, y_hi_scaled, tn_scaled, tau, s, value_vec, discount){
  yhi <- y_hi_scaled - mu * tn_scaled
  ylo <- y_lo_scaled - mu * tn_scaled
  p <- pnorm(yhi, mu, s) - pnorm(ylo, mu, s)
  discount * mu + (pnorm(ylo[1], mu, s) * value_vec[1] +
                     (1 - pnorm(yhi[length(yhi)], mu, s)) * value_vec[length(value_vec)] +
                     sum(p * value_vec))
}

# This version (2018-07) adds an extra block of states which have V calculated directly using only the mu and no
# further learning. The values of these states are used in the calculation of values for
# the other states. The new states are have higher mu values: the original states are within xi
# standard deviations, the new states within xi + EXTRA standard deviations. EXTRA=1 works well.

#' Value calculation for the one-armed bandit with Normal rewards.
#'
#' Assumes `Sigma=mu=0`.
#'
#' @inheritParams nmab_v_args
#' @inheritParams nmab_args
#'
#' @return Difference in value between safe and unknown arms
#'
#' @export
#'
nmab_gi_value <- function(lambda, n, gamma, tau, N, xi, delta){
  EXTRA <- 1 # number of extra xi used for new states
  h <- N + 1
  delta <- delta / sqrt(n) # adjust delta so number of states is constant with n
  mu_range <- seq(0, (xi + EXTRA) * sqrt(1 / n), by=delta)
  lr <- length(mu_range)
  lr2 <- length(seq(0, (xi) * sqrt(1 / n), by=delta))
  value <- matrix(nrow=h, ncol=lr)
  # Value of end states (at stage N)
  value[h, ] <- pmax(mu_range, lambda) * gamma ^ N / (1 - gamma)
  rr <- gamma ^ (0 : N) / (1 - gamma)
  value[, (lr2 : lr)] <- matrix(rr, ncol = 1) %*% matrix(mu_range[lr2 : lr], nrow = 1)
  lo <- mu_range - delta / 2 # lower end of discrete approximation to mu
  hi <- mu_range + delta / 2
  for (j in N : 2){
    t <- j - 1
    tn <- n + tau * (j - 1)
    #the next 3 variables are used for speed-up only
    y_hi_scaled <- hi * (tn + tau) / tau
    y_lo_scaled <- lo * (tn + tau) / tau
    tn_scaled <- tn / tau
    s <- sqrt(1 / tn + 1 / tau) #sd of y
    discount <- gamma ^ t
    safe_reward <- lambda * discount / ( 1 - gamma)
    value_vec <- value[j + 1, ]
    for (i in (lr2 - 1) : 1){
      risky_reward <- nmab_risky_reward(mu_range[i], y_lo_scaled, y_hi_scaled, tn_scaled,
                                        tau, s, value_vec, discount)
      if (risky_reward > safe_reward){
        value[j, i] <- risky_reward
      }else{
        value[j, 1 : i] <- safe_reward
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
  value_vec=value[2, ]
  risky_reward <- nmab_risky_reward(mu_range[1], y_lo_scaled = lo * (n + tau) / tau, y_hi_scaled = hi * (n + tau) / tau,
                                    tn_scaled = n / tau, tau, s, value_vec, discount=1)
  return(risky_reward - lambda / (1 - gamma))
}
