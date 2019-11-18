#' Fits an unbiased logistic regression model 
#' 
#' @param X design matrix
#' @param Y response variable
#' @param estimate_gamma if TRUE, signal strength will be estimated from data
#' 
#' @export 
#' 
reg_log_hd <- function(X, Y, estimate_gamma = FALSE) { 
  
  get_gamma_estimator()
  n_observations <- length(Y)
  n_variables <- ncol(X)
  kappa <- n_variables/n_observations
  solution <-  as.numeric(solve_nonlinear_equations(kappa, gamma))
  alpha_star <- solution[1]
  sigma_star <- solution[2]
  mle_estimator <- get_mle_estimator(Y, X, n_variables)
  c_as <- calculate_c_os(n_observations, n_variables, gamma, 
                         alpha_star, sigma_star)
  
  as_estimator(mle_estimator, c_as)
} 

