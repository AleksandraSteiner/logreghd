library(sigmoid)

calculate_c_os <- function(n_observations, n_variables, gamma, 
                           alfa_star = 1.1678, sigma_star = 3.3466) {
  numerator <- n_observations*gamma*alfa_star
  denominator <-((sigma_star^2)*n_variables + n_observations*gamma*(alfa_star^2))
  numerator/denominator
}

generate_beta_vector <- function(n_variables, n_non_zero, gamma) {
  amplitude <- sqrt(gamma)/(sqrt(n_non_zero))
  c(rep(amplitude, times = n_non_zero), 
    rep(0, times = n_variables - n_non_zero)) 
}

generate_betas <- function(n_observations, kappa, non_zero_proportion, gamma) {
  # TODO: reproducible!
  n_variables <- round(kappa*n_observations)

  beta <- generate_beta_vector(n_variables, 
                               round(non_zero_proportion*n_variables), 
                               gamma)
  
  X <- matrix(rnorm(n_observations*n_variables)/sqrt(n_observations),
              n_observations, n_variables)
  P <- sigmoid(X%*%beta)
  Y <- rbinom(n_observations, 1, P)
  model <- glm(Y ~ X, family = binomial)
  
  list(true = beta,
       mle = model$coef[2:(n_variables+1)])
}

bias <- function(true, estimated) {
  mean(true - estimated)
}
mse <- function(true, estimated) {
  mean((true - estimated)^2)
}
new_estimator <- function(true, estimated, multiplier) {
  mean((true - estimated)*multiplier)
}

simulation <- function(num_iters, n_observations, kappa, non_zero_proportion,
                       gamma, alfa_star = 1.1678) {
  c_os <- calculate_c_os(n_observations, round(kappa*n_observations), 
                         kappa, gamma)
  # TODO: reproducible!
  simulations <- lapply(1:num_iters, function(iteration) {
      betas <- generate_betas(n_observations, kappa, non_zero_proportion, gamma)
      
      list(iteration = iteration,
                 bias_mle = bias(betas[["true"]], betas[["mle"]]),
                 mse_mle = mse(betas[['true']], betas[['mle']]),
                 new_mle = new_estimator(betas[['true']], betas[['mle']], 
                                         betas[["true"]]),
                 bias_ec = bias(betas[['true']]/alfa_star, betas[['mle']]),
                 mse_ec = mse(betas[['true']]/alfa_star, betas[['mle']]),
                 new_ec = new_estimator(betas[['true']]/alfa_star, 
                                        betas[['mle']], 
                                        betas[["true"]]),
                 bias_os = bias(betas[['true']]*c_os, betas[['mle']]),
                 mse_os = mse(betas[['true']]*c_os, betas[['mle']]),
                 new_os = new_estimator(betas[['true']]*c_os, betas[['mle']], 
                                        betas[["true"]]))     
  })
  simulations
}


simulation2 <- function(num_iters, n_observations, kappa, non_zero_proportion,
                       gamma, alfa_star = 1.1678) {
  c_os <- calculate_c_os(n_observations, round(kappa*n_observations), 
                         kappa, gamma)
  # TODO: reproducible!
  simulations <- lapply(1:num_iters, function(iteration) {
    n_variables <- round(kappa*n_observations)
    
    beta <- generate_beta_vector(n_variables, 
                                 round(non_zero_proportion*n_variables), 
                                 gamma)
    
    X <- matrix(rnorm(n_observations*n_variables)/sqrt(n_observations),
                n_observations, n_variables)
    P <- sigmoid(X%*%beta)
    Y <- rbinom(n_observations, 1, P)
    model <- glm(Y ~ X, family = binomial)
    
    betas <- list(true = beta,
                  mle = model$coef[2:(n_variables+1)])
    
    
    list(iteration = iteration,
         bias_mle = bias(betas[["true"]], betas[["mle"]]),
         mse_mle = mse(betas[['true']], betas[['mle']]),
         new_mle = new_estimator(betas[['true']], betas[['mle']], 
                                 betas[["true"]]),
         bias_ec = bias(betas[['true']]/alfa_star, betas[['mle']]),
         mse_ec = mse(betas[['true']]/alfa_star, betas[['mle']]),
         new_ec = new_estimator(betas[['true']]/alfa_star, 
                                betas[['mle']], 
                                betas[["true"]]),
         bias_os = bias(betas[['true']]*c_os, betas[['mle']]),
         mse_os = mse(betas[['true']]*c_os, betas[['mle']]),
         new_os = new_estimator(betas[['true']]*c_os, betas[['mle']], 
                                betas[["true"]]))     
  })
  simulations
}
