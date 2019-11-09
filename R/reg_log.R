#' Fits an unbiased logistic regression model 
#' 
#' @param X design matrix
#' @param Y response variable
#' @param estimate_gamma if TRUE, signal strength will be estimated from data
#' 
#' @export 
#' 
reg_log_hd <- function(X, Y, estimate_gamma = FALSE) { 
  
  if(estimate_gamma == FALSE) {
    gamma <- 5
  }
  else {
    stop('not implemented yet')
  }
  
  n_observations <- length(Y)
  n_variables <- length(X[1, ])
  kappa <- n_variables/n_observations
  alpha_star <- solve_nonlinear_equations(kappa, gamma)$alpha
  sigma_star <- solve_nonlinear_equations(kappa, gamma)$sigma
  
  logistic_model <- glm(Y ~ X, family = binomial)
  mle_estimator <- logistic_model$coefficients[2 : (n_variables+1)]
  c_as <- calculate_c_os(n_observations, n_variables, gamma, 
                         alpha_star, sigma_star)
  
  as_estimator(mle_estimator, c_as)
} 