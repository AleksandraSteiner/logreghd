calculate_c_os <- function(n_observations, n_variables, gamma, 
                           alfa_star = 1.1678, sigma_star = 3.3466) {
  numerator <- n_observations*gamma*alfa_star
  denominator <-((sigma_star^2)*n_variables + n_observations*gamma*(alfa_star^2))
  numerator/denominator
}

as_estimator <- function(mle_estimator, c_as) {
  mle_estimator * c_as
}

get_kappa_grid <- function(method, n) {
  if(method == 'random') {
    runif(n, 0, 0.5)
  } else {
    stop('not implemented yet')
  }
}

mle_existance_bootstrap <- function(X, kappa, n_rep) {  #n_rep > 0
  lapply(1:n_rep, function(iteration) {
    mle_existance(X, round(ncol(X)/kappa))
  })
}

mle_existance <- function(X, n_obs) {
  bootstrap_sample <- sample(1:nrow(X), n_obs, replace = TRUE)
  bootstrap_data <- X[bootstrap_sample, ]
  

  }

estimate_gamma <- function(X, Y, n_kappa = 50, kappa_method = 'random', n_rep = 10) {
  kappa_grid <- get_kappa_grid(kappa_method, n_kappa) #TODO: add minimum kappa as parameter
  
  lapply(kappa_grid, function(kappa) {
    
  })
}