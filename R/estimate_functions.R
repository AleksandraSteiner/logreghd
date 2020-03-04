#' @importFrom stats runif
get_kappa_grid <- function(method, n_kappa, kappa) {
  if(method == 'random') {
    runif(n_kappa, kappa, 0.5)
  } else {
    seq(kappa, 0.5, length.out = n_kappa)
  }
}

#' @importFrom lpSolve lp
#checks if MLE exists
is_feasible <- function(X, Y) { 
  solution_object <- lp("min", 0, X, 
                        ifelse(Y == 0, "<=", ">="), 
                        ifelse(Y == 0, -1, 1)) 
  #we do not put zero because of weak inequality
  solution_object[28] == 0 
  #if true, the problem is feasible = MLE does not exist
}

#calculates how often MLE doesn't exist for i-th kappa, 
#so it calculates pi for one kappa
estimate_pi_i <- function(X, Y, n_rep, kappa_i) { 
  n_variables <- ncol(X)
  n_observations <- nrow(X)
  feasibile <- sapply(1:n_rep, function(i) {
    ind <- sample(1:n_observations, round(n_variables/kappa_i)) 
    is_feasible(X[ind, ], Y[ind])
  })
  sum(feasibile)/n_rep
}

#caluculates pi_hats for kappa grid
estimate_pi <- function(X, Y, n_rep, kappa_grid) {
  n_variables <- ncol(X)
  n_observations <- nrow(X)
  sapply(kappa_grid, function(kappa_i) {
    estimate_pi_i(X, Y, n_rep, kappa_i)
  })
}

#finds a pair of boundary kappas to find then such kappa for which pi = 0.5 
find_boundary <- function(kappa_grid, pi_grid) {
  j <- which(pi_grid >= 0.5)[1]
  if(j > 1) {
    list(kappa_i = kappa_grid[j - 1], 
         kappa_j = kappa_grid[j], 
         pi_hat_kappa_i = pi_grid[j - 1], 
         pi_hat_kappa_j = pi_grid[j])
  } else if (j == 1) {
    else if(pi_grid[1] == 0.5) {
      list(kappa_i = kappa_grid[j], 
           kappa_j = kappa_grid[j], 
           pi_hat_kappa_i = pi_grid[j], 
           pi_hat_kappa_j = pi_grid[j] + 1)
    }
    else if(pi_grid[1] > 0.5) {
      #Is it possible and should be kappa <- kappa_grid[1] found ?
    }
  } else if (j == integer(0)) {
    stop('no kappa estimate exists, try with a finer grid')  } #is it a correct message ?
}

#' @importFrom stats approx
estimate_kappa <- function(X, Y, n_rep = 10, 
                           n_kappa = 50, 
                           kappa_method = 'random') {
  n_variables <- ncol(X)
  n_observations <- nrow(X)
  kappa_grid <- get_kappa_grid(kappa_method, n_kappa, 
                               n_variables/n_observations) 
  pi_hat_kappa <- estimate_pi(n_rep, n_observations, 
                                    n_variables, kappa_grid)
  boundary <- find_boundary(kappa_grid, pi_hat_kappa)
  kappa_j_1 <- boundary[1]
  kappa_j <- boundary[2]
  pi_j_1 <- boundary[3]
  pi_j <- boundary[4]
  kappa_j_1 + (0.5 - pi_j_1) * (kappa_j - kappa_j_1) / (pi_j - pi_j_1) 
  #interpolation solution for kappa_hat (x)
}

#' @importFrom R.matlab readMat
estimate_gamma <- function(kappa_hat) {
  kappa <- kappa_gamma_data[,1] == round(kappa_hat, digits = 3) #round/floor ?
  kappa_gamma_data[kappa, 2]
}

#two-dimensional density
calculate_Q_density <- function(q1, q2, alpha, gamma, kappa, sigma) {
  exp((-q1^2 * ( (alpha * gamma)^2 + kappa * sigma^2)  
       - 2 * q1 * q2 * alpha * gamma^2
       - q2^2 * gamma^2) / 2 * kappa * sigma^2) 
  / (2 * pi * sigma * sqrt(kappa))
}

prox_function <- function(z) {
  
}

#functions from system of equations
h_function <- function(q1, q2) {
  (2 * lambda^2 * exp(q1) * exp(2 * prox(q2))) /
    ((1 + exp(q1)) * (1 + exp(prox(q2)))^2)
}

g_function <- function(q1, q2) {
  (lambda * q1 * exp(q1) * exp(prox(q2))) /
    ((1 + exp(q1)) * (1 + exp(prox(q2))))
}

j_function <- function(q1, q2) {
  (2 * exp(q1) * (1 + exp(prox(q2)))^2) /
    ((1 + exp(q1)) * ((1 + exp(prox(q2)))^2 + lambda * exp(prox(q2))))
}

calculate_values <- function(alpha, sigma, lambda) {
  y <- numeric(3)
  y[1] <- (1/kappa^2) * adaptIntegrate(h_function, 
                               lowerLimit = c(-infinity, -infinity), 
                               upperLimit = c(+infinity, +infinity)) - sigma^2
  y[2] <- adaptIntegrate(g_function, 
                         lowerLimit = c(-infinity, -infinity), 
                         upperLimit = c(+infinity, +infinity))
  y[3] <- adaptIntegrate(h_function, 
                         lowerLimit = c(-infinity, -infinity), 
                         upperLimit = c(+infinity, +infinity)) + kappa - 1
  y
}

#' @importFrom nleqslv nleqslv
solve_equations <- function(kappa, gamma) {
  xstart <- c() #what should be the starting point here ?
  solution <- nleqslv(xstart, calculate_values)$x #control options
  alpha <- solution[1]
  sigma <- solution[2]
  lambda <- solution[3]
  list(alpha = alpha, sigma = sigma, lambda = lambda)
}

get_gamma_estimator <- function(X, Y, estimate_gamma) {
  if(estimate_gamma == FALSE) {
    gamma <- 5
  }
  else {
    stop('not implemented yet')
  }
  gamma
}

#' @importFrom stats glm binomial coef
get_mle_estimator <- function(Y, X, n_variables) {
  logistic_model <- glm(Y ~ X, family = binomial)
  coef(logistic_model)[2 : (n_variables + 1)]
}

calculate_c_os <- function(n_observations, n_variables, gamma, 
                           alfa_star = 1.1678, sigma_star = 3.3466) {
  numerator <- n_observations*gamma*alfa_star
  denominator <-((sigma_star^2)*n_variables + n_observations*gamma*(alfa_star^2))
  numerator/denominator
}

get_as_estimator <- function(mle_estimator, c_as) {
  mle_estimator * c_as
}

get_ec_estimator <- function(mle_estimator, alpha_star) {
  mle_estimator / alpha_star
}