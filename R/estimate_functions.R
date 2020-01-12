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

#checks if MLE exists
is_feasible <- function(X, Y) { 
  solution_object <- lp("min", 0, X, 
                        ifelse(Y == 0, "<=", ">="), 
                        ifelse(Y == 0, -1, 1)) 
                        #we do not put zero because of weak inequality
  solution_object[28] == 0 
  #if true, the problem is feasible, so MLE does not exist
}

#caclulates how often MLE doesn't exist for i-th kappa
estimate_pi_i <- function(X, Y, n_rep, kappa_i) { 
  n_variables <- ncol(X)
  n_observations <- nrow(X)
  feasibile <- sapply(1:n_rep, function(i) {
    ind <- sample(1:n_observations, round(n_variables/kappa_i)) 
    is_feasible(X[ind, ], Y[ind])
  })
  sum(feasibile)/n_rep
  print("hello world")
}

#caluculates pi_hats for kappa grid
estimate_pi_kappagrid <- function(X, Y, n_rep, kappa_grid) {
  n_variables <- ncol(X)
  n_observations <- nrow(X)
  sapply(kappa_grid, function(kappa_i) {
    estimate_pi_i(X, Y, n_rep, kappa_i)
  })
}

find_boundary <- function(kappa_grid, pi_kappas) {
  j <- which(pi_kappas >= 0.5)[1]
  if(j > 1) {
    list(kappa_i = kappa_grid[j - 1], 
         kappa_j = kappa_grid[j], 
         pi_hat_kappa_i = pi_kappas[j - 1], 
         pi_hat_kappa_j = pi_kappas[j])
  } else if (j == 1) {
    #use new kappa < kappa_1 and corresponding pi ?
  } else if (j == integer(0)) {
    #return proper information/suggestion if there is not such kappa ?
  }
}

estimate_kappa <- function(X, Y, n_rep = 10, 
                           n_kappa = 50, 
                           kappa_method = 'random') {
  n_variables <- ncol(X)
  n_observations <- nrow(X)
  kappa_grid <- get_kappa_grid(kappa_method, n_kappa, 
                               n_variables/n_observations) 
  pi_hat_kappa <- estimate_pi_kappa(n_rep, n_observations, 
                                    n_variables, kappa_grid)
  boundary <- find_boundary(kappa_grid, pi_hat_kappa)
  kappa_j_1 <- boundary[1]
  kappa_j <- boundary[2]
  pi_j_1 <- boundary[3]
  pi_j <- boundary[4]
  kappa_j_1 + (0.5 - pi_j_1) * (kappa_j - kappa_j_1) 
  / (pi_j - pi_j_1) 
  #interpolation solution for kappa_hat (x)
}

estimate_gamma <- function(kappa_hat) {
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
  mle_estimator <- logistic_model[[1]][2 : (n_variables + 1)]
  return(mle_estimator)
}