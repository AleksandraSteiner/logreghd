calculate_c_os <- function(n_observations, n_variables, gamma, 
                           alfa_star = 1.1678, sigma_star = 3.3466) {
  numerator <- n_observations*gamma*alfa_star
  denominator <-((sigma_star^2)*n_variables + n_observations*gamma*(alfa_star^2))
  numerator/denominator
}

as_estimator <- function(mle_estimator, c_as) {
  mle_estimator * c_as
}

get_kappa_grid <- function(method, n_kappa, kappa) {
  if(method == 'random') {
    runif(n_kappa, kappa, 0.5)
  } else {
    seq(kappa, 0.5, length.out = n_kappa)
  }
}

is_feasible <- function(X, Y) {
  solution_object <- lp("min", 0, X, 
                        ifelse(Y == 0, "<=", ">="), 
                        ifelse(Y == 0, -1, 1)) #weak inequality
  solution_object[28] == 0
}

estimate_pi <- function(n_rep, n_observations, n_variables, kappa) {
  feasibile <- sapply(1:n_rep, function(i) {
    ind <- sample(1:n_observations, round(n_variables/kappa))
    is_feasible(X[ind, ], Y[ind])
  })
  sum(feasibile)/n_rep
}

estimate_pi_kappa <- function(X, Y, n_rep, n_observations, n_variables, kappa_grid) {
  lapply(kappa_grid, function(kappa) {
    estimate_pi(n_rep, n_observations, n_variables, kappa)
  })
}

find_brink <- function(kappa_grid, pi_hat_kappa) { #TODO what if there is no such kappa
  j <- (pi_hat_kappa >= 0.5)[1]
  list(kappa_i = kappa_grid[j-1], kappa_j = kappa_grid[j], 
       pi_hat_kappa_i = pi_hat_kappa[j-1], pi_hat_kappa_j = pi_hat_kappa[j])
}

estimate_kappa <- function(X, Y, n_rep = 10, n_kappa = 50, kappa_method = 'random') {
  n_variables <- ncol(X)
  n_observations <- nrow(X)
  kappa_grid <- get_kappa_grid(kappa_method, n_kappa, n_variables/n_observations) 
  pi_hat_kappa <- estimate_pi_kappa(n_rep, n_observations, n_variables, kappa_grid)
  brink_points <- find_brink(kappa_grid, pi_hat_kappa)
  interpolation <- approx(c(brink_points[1], brink_points[2]), 
                          c(brink_points[3], brink_points[4]))
  ind <- interpolation[['y']] == 0.5
  x <- interpolation[['x']]
  x[ind]
}

estimate_gamma <- function(X, Y, n_kappa = 50, kappa_method = 'random', n_rep = 10) {
  #TODO: set gamma_hat = g_MLE(kappa_hat)
}

solve_equations <- function(kappa, gamma) {
  alpha <- 1
  sigma <- 1
  lambda <- 1
  list(alpha = alpha, sigma = sigma, lambda = lambda)
}

get_gamma_estimator <- function(X, Y, estimate_gamma) {
  if(estimate_gamma == FALSE) {
    gamma <- 5
  }
  else {
    stop('not implemented yet')
  }
  return(gamma)
}

get_mle_estimator <- function(Y, X, n_variables) {
  logistic_model <- glm(Y ~ X, family = binomial)
  mle_estimator <- logistic_model[1][2 : (n_variables+1)]
  return(mle_estimator)
}