#' Calculate Gittins indices for multiple arms (Bernoulli rewards)
#'
#' The states are a triangular matrix with:
#' Sigma = Sigma_start, Sigma_start + 1, ..., Sigma_start + num_actions - 1,
#' n = n_start, n_start + 1, ..., n_start + num_actions - 1,
#' and Sigma <= Sigma_start + n - n_start.
#'
#' @param Sigma_start=1 lowest value of Sigma for the arms
#' @param n_start=2 lowest value of n for the arms
#' @param gamma numeric in (0, 1]; discount factor
#' @param N integer>0; time horizon used
#' @param num_actions determinines the number of states GI are calculated for
#' @param tol absolute accuracy required
#'
#' @return A triangular matrix of GI values
#'
#' @export
#'
bmab_gi_multiple <- function(Sigma_start=1, n_start=2, gamma, N, num_actions, tol){
  GI <- bmab_gi_multiple_ab(alpha_start=Sigma_start, beta_start=n_start - Sigma_start,
                            gamma, N, num_actions, tol)
  for(i in 2 : num_actions){
    GI[, i] <- c(rep(NA, i - 1), GI[1 : (num_actions - i + 1), i])
  }
  GI
}

#' Calculate Gittins indices for multiple arms (Bernoulli rewards)
#'
#' This version uses alpha, beta parameterisation.
#' The states are a triangular matrix with:
#' alpha=alpha_start, alpha_start + 1, ..., alpha_start + num_actions - 1,
#' beta=beta_start, beta_start + 1, ..., beta_start + num_actions - 1,
#' and alpha + beta <= num_actions + alpha_start + beta_start - 1.
#'
#' @param alpha_start=1 lowest value of alpha for the arms
#' @param beta_start=1 lowest value of beta for the arms
#' @param gamma discount factor
#' @param N time horizon used
#' @param num_actions determinines the number of states GI are calculated for
#' @param tol absolute accuracy required
#'
#' @return A triangular matrix of GI values
#'
#' @export
#'
bmab_gi_multiple_ab <- function(alpha_start=1, beta_start=1, gamma, N, num_actions, tol){
  GI <- matrix(rep(NA, num_actions * num_actions), nrow=num_actions, ncol=num_actions)
  alpha_range <- alpha_start : (alpha_start + num_actions - 1)
  beta_range <- beta_start : (beta_start + num_actions - 1)
  mu <- alpha_start / (alpha_start + beta_range)
  if (gamma==1){
    lb_vec <- mu
  }else{
    lb_vec <- bmab_kgi(alpha_start, alpha_start + beta_range, gamma)
  }
  cat("Calculating GI values for", as.integer(0.5 * num_actions * (num_actions + 1)), "states\n")
  pb <- txtProgressBar(min = 0, max = num_actions, style = 3)
  for (a in 1 : num_actions){
    ub <- 1
    for (b in 1 : (num_actions - a + 1)){
      GI[b, a] <- bmab_gi_ab(alpha_range[a], beta_range[b], gamma, tol, N, lb=lb_vec[b], ub=ub)
      ub <- GI[b, a]
    }
    lb_vec <- GI[, a]
    setTxtProgressBar(pb, a)
  }
  close(pb)
  GI
}

#' Calculate Gittins indices for a single arm (Bernoulli rewards)
#'
#' The initial interval for calibration are as follows:
#' For lower bound, use lb if supplied else use KGI if kgi=T or Sigma/n otherwise.
#' For upper bound, use ub if supplied else use GI+ if giplus=T or 1 otherwise.
#'
#' @param Sigma value of Sigma for the arm
#' @param n value of n for the arm
#' @param gamma discount factor
#' @param tol absolute accuracy required
#' @param N time horizon used
#' @param lb=NA optional lower bound for GI
#' @param ub=NA optional upper bound for GI
#' @param kgi=F optional boolean indicates whether to use KGI for lower bound (only if lb=NA)
#' @param giplus=F optional boolean indicates whether to use GI+ for upper bound (only if ub=NA)
#'
#' @return A vector of GI values
#'
#' @export
#'
bmab_gi <- function(Sigma, n, gamma, tol, N, lb=NA, ub=NA, kgi=F, giplus=F){
  if (is.na(lb)){
    if (kgi){
      lb <- bmab_kgi(Sigma, n, gamma)
    }else{
      lb <- Sigma / n
    }
  }
  if (is.na(ub)){
    if (giplus){
      ub <- bmab_giplus(Sigma, n, gamma, tol, upper=T)
    }else{
      ub <- 1
    }
  }
  mean(calibrate_arm(bmab_gi_value, lb, ub, tol, Sigma, n, gamma, N))
}

#' Calculate Gittins indices for a single arm (Bernoulli rewards)
#'
#' This version uses alpha, beta parameterisation.
#' The initial interval for calibration are as follows:
#' For lower bound, use lb if supplied else use KGI if kgi=T or Sigma/n otherwise.
#' For upper bound, use ub if supplied else use GI+ if giplus=T or 1 otherwise.
#'
#' @param alpha value of alpha for the arm
#' @param beta value of beta for the arm
#' @param gamma discount factor
#' @param tol absolute accuracy required
#' @param N time horizon used
#' @param lb=NA optional lower bound for GI
#' @param ub=NA optional upper bound for GI
#' @param kgi=F optional boolean indicates whether to use KGI for lower bound (only if lb=NA)
#' @param giplus=F optional boolean indicates whether to use GI+ for upper bound (only if ub=NA)
#'
#' @return A vector of GI values
#'
#' @export
#'
bmab_gi_ab <- function(alpha, beta, gamma, tol, N, lb=NA, ub=NA, kgi=F, giplus=F){
  bmab_gi(Sigma = alpha, n = alpha + beta, gamma, tol, N, lb, ub, kgi, giplus)
}

#' Calculate the GI+ index for a single arm (Bernoulli rewards)
#'
#' Upper bound for GI.
#'
#' @param Sigma value of Sigma for the arm
#' @param n value of n for the arm
#' @param gamma discount factor
#' @param tol absolute accuracy required
#' @param upper=F if TRUE, the upper end of the interval is returned, otherwise the midpoint
#'
#' @return A vector of GI+ values
#'
#' @export
#'
bmab_giplus <- function(Sigma, n, gamma, tol, upper=F){
  interval <- calibrate_arm(bmab_giplus_value, lb=Sigma / n, ub=1, tol, Sigma, n, gamma)
  if (upper){
    return(interval[2])
  }
  mean(interval)
}

#' Calculate the knowledge gradient index for a single arm (Bernoulli rewards)
#'
#' Exact closed form calculation.
#'
#' @param Sigma value of Sigma for the arm
#' @param n value of n for the arm
#' @param gamma discount factor
#'
#' @return A vector of KGI values
#'
#' @export
#'
bmab_kgi <- function(Sigma, n, gamma){
  mu <- Sigma / n
  H <- gamma / (1 - gamma)
  (mu + H * mu * (Sigma + 1) / (n + 1)) / (1 + H * mu)
}

#' Value of one-armed bandit using GI+ (Bernoulli rewards)
#'
#' @param lambda reward from the known arm
#' @param Sigma value of Sigma for the unknown arm
#' @param n value of n for the unknown arm
#' @param gamma discount factor
#'
#' @return Difference in value between safe and unknown arms
#'
#' @export
#'
bmab_giplus_value <- function(lambda, Sigma, n, gamma){
  mu <- Sigma / n
  mu_success <- (Sigma + 1) / (n + 1)
  H <- gamma / (1 - gamma)
  continue <- function(x) {dbeta(x, Sigma + 1, n - Sigma) * x}
  value_success <- H * (integrate(continue, lambda, 1)[[1]] + lambda * pbeta(lambda, Sigma + 1, n - Sigma))
  value_fail <- H * lambda
  mu + mu * value_success + (1 - mu) * value_fail - lambda / (1 - gamma)
}

#' Value of one-armed bandit using KGI (Bernoulli rewards)
#'
#' @param lambda reward from the known arm
#' @param Sigma value of Sigma for the unknown arm
#' @param n value of n for the unknown arm
#' @param gamma discount factor
#' @param N time horizon used in the dynamic programme
#'
#' @return Difference in value between safe and unknown arms
#'
#' @export
#'
bmab_gi_value <- function(lambda, Sigma, n, gamma, N){
  h <- N + 1
  n_vec <- n : (n + N)
  s_vec <- Sigma : (Sigma + N)
  mu <- outer(s_vec, n_vec, "/")
  value_mat <- matrix(nrow=h, ncol=h)
  # Values of end states
  if (gamma==1){
    value_mat[, h] <- pmax(mu[, h], lambda)
    safe_reward <- lambda * (N + 2 - (1 : N))
  }else{
    value_mat[, h] <- pmax(mu[, h], lambda) * gamma ^ N / (1 - gamma)
    safe_reward <- lambda * gamma ^ ((1 : N) - 1) / (1 - gamma)
  }
  # Run DP to get values of other states
  for (i in N : 1){
    j <- i + 1
    risky_reward <- mu[1 : i, i] * (gamma ^ (i - 1) + value_mat[2 : j, j]) +
      (1 - mu[1 : i, i]) * value_mat[1 : i, j]
    value_mat[1 : i, i] <- pmax(risky_reward, safe_reward[i])
  }
  return(value_mat[1, 1] - lambda / (1 - gamma))
}
