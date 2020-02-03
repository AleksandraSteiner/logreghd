#' @title Fits MSE-minimizing logistic regression model 
#' 
#' @description give a deeper explanation
#' 
#' @param X design matrix
#' @param Y response variable
#' @param estimate_gamma if TRUE, signal strength will be estimated from data
#' 
#' @seealso \code{\link[stats]{glm}}
#' 
#' @export 
#' 
#' 
reg_log_hd <- function(X, Y, estimate_gamma = FALSE) { 
  
  n_observations <- nrow(Y)
  n_variables <- ncol(X)
  gamma <- get_gamma_estimator(X, Y, estimate_gamma)
  solution <-  solve_equations(n_variables/n_observations, gamma)
  
  as_estimator(get_mle_estimator(Y, X, n_variables), 
               calculate_c_os(n_observations, n_variables, gamma, 
                              solution[[1]], solution[[2]]))
}