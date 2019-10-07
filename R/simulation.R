# data: 100x10
# petla: 1000

alfa_star <- 1.1678
sigma_star <- 3.3466
ile <- 1000

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

cdiff <- function(true_beta, estimated_beta, transformation = identity,
                  mle_multiplier = 1) {
  mean(transformation(estimated_beta*mle_multiplier - true_beta))
}

simulation <- function(num_iters, n_observations, kappa, non_zero_proportion,
                       gamma, alfa_star = 1.1678) {
  c_os <- calculate_c_os(n_observations, round(kappa*n_observations), 
                         kappa, gamma)
  # TODO: reproducible!
  simulations <- lapply(1:num_iters, function(iteration) {
      betas <- generate_betas(n_observations, kappa, non_zero_proportion, gamma)
      
      data.frame(iteration = iteration,
                 bias_mle = cdiff(betas[["true"]], betas[["mle"]]),
                 mse_mle = cdiff(betas[['true']], betas[['mle']], function(x) x^2),
                 new_mle = cdiff(betas[['true']], betas[['mle']], 
                                 function(x) x*betas[["true"]]),
                 bias_ec = cdiff(betas[['true']], betas[['mle']], 
                                 mle_multiplier = 1/alfa_star),
                 mse_ec = cdiff(betas[['true']], betas[['mle']], 
                                function(x) x^2,
                                mle_multiplier = 1/alfa_star),
                 new_ec = cdiff(betas[['true']], betas[['mle']], 
                                function(x) x*betas[["true"]],
                                mle_multiplier = 1/alfa_star),
                 bias_os = cdiff(betas[['true']], betas[['mle']], 
                                 mle_multiplier = c_os),
                 mse_os = cdiff(betas[['true']], betas[['mle']], 
                                function(x) x^2,
                                mle_multiplier = c_os),
                 new_os = cdiff(betas[['true']], betas[['mle']], 
                                function(x) x*betas[["true"]],
                                mle_multiplier = c_os))     
  })
  # TODO: remove_dplyr
  dplyr::bind_rows(simulations)
}
# 
# lapply(
#   c(idenity, function(x) x^2) {
#     lapply(multipliers, {
#       ...
#     } )
#   }
# ) %>%
#   unlist()

proba <- simulation(num_iters = 100, n_observations = 1000, kappa = 0.2, 
                    non_zero_proportion = 0.5, gamma = 5)
proba
