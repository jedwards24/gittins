library(bench)

tt <- bmab_gi_multiple_ab(1, 1, gamma = 0.5, N = 100, num_actions = 4, tol = 1e-5)
tt
as.vector(tt)

bmab_gi_multiple_ab3 <- function(alpha_start = 1, beta_start = 1, gamma, N, num_actions, tol){
  n_states <- as.integer(round(0.5 * num_actions * (num_actions + 1), 0))
  alpha_vec <- beta_vec <- integer(0)
  GI <- matrix(rep(NA, num_actions * num_actions), nrow = num_actions, ncol = num_actions)
  alpha_range <- alpha_start:(alpha_start + num_actions - 1)
  beta_range <- beta_start:(beta_start + num_actions - 1)
  mu <- alpha_start / (alpha_start + beta_range)
  lb_vec <- bmab_kgi(alpha_start, alpha_start + beta_range, gamma)
  cat("Calculating GI values for", as.integer(0.5 * num_actions * (num_actions + 1)), "states\n")
  pb <- txtProgressBar(min = 0, max = num_actions, style = 3)
  for (a in 1:num_actions){
    ub <- 1
    betas <- beta_range[1:(num_actions - a + 1)]
    alpha_vec <- c(alpha_vec, rep(a, length(betas)))
    beta_vec <- c(beta_vec, betas)
    for (b in 1:(num_actions - a + 1)){
      GI[b, a] <- bmab_gi_ab(alpha_range[a], beta_range[b], gamma, tol, N, lb = lb_vec[b], ub = ub)
      ub <- GI[b, a]
    }
    lb_vec <- GI[, a]
    setTxtProgressBar(pb, a)
  }
  close(pb)
  data.frame(alpha = alpha_vec, beta = beta_vec, gi = as.numeric(na.omit(as.vector(GI))))
}

bmab_gi_multiple_ab4 <- function(alpha_start = 1, beta_start = 1, gamma, N, num_actions, tol){
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
      gi_mat[b, a] <- bmab_gi_ab(alpha_range[a], beta_range[b], gamma, tol, N, lb = lb_vec[b], ub = ub)
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
  df <- tibble::tibble(alpha = alpha_vec, beta = beta_vec, gi = as.numeric(na.omit(as.vector(gi_mat))),
                   Sigma = alpha_vec, n = n_vec, stage = n_vec - min(n_vec))
  attr(df, "params") <- list(alpha_start = alpha_start, beta_start = beta_start,
                             Sigma_start = alpha_start, n_start = alpha_start + beta_start,
                             gamma = gamma, N = N, num_actions = num_actions, tol = tol)
  attr(df, "gi_matrix") <- gi_mat
  attr(df, "gi_matrix_ns") <- gi_mat_ns
  df
}

df <- bmab_gi_multiple_ab3(1, 1, gamma = 0.5, N = 100, num_actions = 4, tol = 1e-5)
df2 <- bmab_gi_multiple_ab4(1, 1, gamma = 0.5, N = 100, num_actions = 4, tol = 1e-5)
identical(df, df2)
df2
attributes(df2)
f <- function(x = 1, y = 2) {
  as.list(...)
}
f()
df$gi
as.vector(df$gi)

as.vector(na.omit(tt))
as.numeric(na.omit(as.vector(tt)))

bench::mark(f3 = bmab_gi_multiple_ab3(1, 1, gamma = 0.5, N = 100, num_actions = 4, tol = 1e-5),
            f4 = bmab_gi_multiple_ab4(1, 1, gamma = 0.5, N = 100, num_actions = 4, tol = 1e-5))
tt
system_time(bmab_gi_multiple_ab3(1, 1, gamma = 0.5, N = 100, num_actions = 20, tol = 1e-5))
system_time(bmab_gi_multiple_ab4(1, 1, gamma = 0.5, N = 100, num_actions = 20, tol = 1e-5))

test <- function(alpha_start = 1, beta_start = 1, gamma, N, num_actions, tol){
  alpha_vec <- beta_vec <- integer(0)
  alpha_range <- alpha_start:(alpha_start + num_actions - 1)
  beta_range <- beta_start:(beta_start + num_actions - 1)
  for (a in 1:num_actions){
    betas <- beta_range[1:(num_actions - a + 1)]
    alpha_vec <- c(alpha_vec, rep(a, length(betas)))
    beta_vec <- c(beta_vec, betas)
  }
  data.frame(alpha = alpha_vec, beta = beta_vec)
}
test(1, 1, gamma = 0.5, N = 100, num_actions = 4, tol = 1e-5)
tt
