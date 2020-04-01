#' @importFrom stats runif
get_kappa_grid <- function(method, n_kappa, kappa) {
  if(method == "random") {
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
  } 
  else if (j == 1) {
    if (pi_grid[1] == 0.5) {
      list(kappa_i = kappa_grid[j], 
           kappa_j = kappa_grid[j], 
           pi_hat_kappa_i = pi_grid[j], 
           pi_hat_kappa_j = pi_grid[j] + 1)
    }
    else if (pi_grid[1] > 0.5) { #down interpolation
      list(kappa_i = kappa_grid[j], 
           kappa_j = kappa_grid[j + 1], 
           pi_hat_kappa_i = pi_grid[j], 
           pi_hat_kappa_j = pi_grid[j + 1])
    }
  }
  else if (j == integer(0)) {
    stop("No kappa estimate exists, try with a finer grid")  } 
}

estimate_kappa <- function(X, Y, n_rep = 10, 
                           n_kappa = 50, 
                           kappa_method = "random") {
  n_variables <- ncol(X)
  n_observations <- nrow(X)
  kappa_grid <- get_kappa_grid(kappa_method, n_kappa, 
                               n_variables / n_observations) 
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

estimate_gamma <- function(kappa_hat) {
  kappa <- kappa_gamma_data[, 1] == round(kappa_hat, digits = 3) #round/floor ?
  kappa_gamma_data[kappa, 2]
}

#two-dimensional density
calculate_Q_density <- function(alpha, sigma, kappa, gamma) {
  density <- function(q1, q2) {
  exp((-q1 ^ 2 * ( (alpha * gamma) ^ 2 + kappa * sigma ^ 2)  
       - 2 * q1 * q2 * alpha * gamma ^ 2 - q2 ^ 2 * gamma ^ 2) / 2 * kappa * sigma ^ 2) /
    (2 * pi * sigma * sqrt(kappa))
  }
}

prox <- function(z, lambda) {
  der <- function(t, z1 = z, lambda1 = lambda) {
    lambda1 * exp(t) / (1 + exp(t)) + t - z1
  }
  derivative <- function(t) {
    xstart <- 1
    solution <- nleqslv(xstart, der)$x 
    solution[1]
  }
  derivative()
}

integral_1 <- function(alpha, sigma, lambda, kappa, gamma) {
  function(q) {
    calculate_Q_density(q[1], q[2], alpha, gamma, kappa, sigma) *
      ((2 * lambda ^ 2 * exp(q[1]) * exp(2 * prox(q[2], lambda))) /
         ((1 + exp(q[1])) * (1 + exp(prox(q[2], lambda))) ^ 2 * kappa ^ 2)) - sigma ^ 2
  }
}

integral_2 <- function(alpha, sigma, lambda, kappa, gamma) {
  function(q) {
    calculate_Q_density(q[1], q[2], alpha, gamma, kappa, sigma) *
      (lambda * q[1] * exp(q[1]) * exp(prox(q[2], lambda))) /
      ((1 + exp(q[1])) * (1 + exp(prox(q[2], lambda))))
  }
}

integral_3 <- function(alpha, sigma, lambda, kappa, gamma) {
  function(q) {
    calculate_Q_density(q[1], q[2], alpha, gamma, kappa, sigma) *
      (2 * exp(q[1]) * (1 + exp(prox(q[2], lambda))) ^ 2) /
      ((1 + exp(q[1])) * ((1 + exp(prox(q[2], lambda))) ^ 2 
                          + lambda * exp(prox(q[2], lambda)))) + kappa - 1
  }
}

#' @importFrom cubature adaptIntegrate
calculate_integral <- function(alpha, sigma, lambda, kappa, gamma, integral) {
  adaptIntegrate(integral, c(-Inf, -Inf), c(+Inf, +Inf))[[1]]
}

system_values <- function(kappa, gamma) {
  expected_values <- function(alpha, sigma, lambda) {
    c(calculate_integral(alpha, sigma, lambda, kappa, gamma, integral_1),
      calculate_integral(alpha, sigma, lambda, kappa, gamma, integral_2),
      calculate_integral(alpha, sigma, lambda, kappa, gamma, integral_3))
  }
  expected_values()
}

#' @importFrom nleqslv nleqslv
solve_equations <- function(kappa, gamma) {
  xstart <- c(1.1678, 3.3466, 0.9605) #these are solutions for kappa = 0.1 and gamma = sqrt(5)
  system_values <- system_values(kappa, gamma)
  solution <- nleqslv(xstart, system_values)$x 
  names(solution) <- c("alpha", "sigma", "lambda")
  solution
}

get_gamma_estimator <- function(X, Y, estimate_gamma) {
  if(estimate_gamma == FALSE) {
    gamma <- sqrt(5)
  }
  else {
    stop("not implemented yet")
  }
  gamma
}

#' @importFrom stats glm binomial coef
get_mle_estimator <- function(Y, X, n_variables) {
  logistic_model <- glm(Y ~ X, family = binomial)
  coef(logistic_model)[2 : (n_variables + 1)]
}

calculate_c_os <- function(n_observations, n_variables, gamma, 
                           alpha_star = 1.1678, sigma_star = 3.3466) {
  numerator <- n_observations * gamma * alpha_star
  denominator <- (sigma_star ^ 2) * n_variables + 
    n_observations * gamma * (alpha_star ^ 2)
  numerator / denominator
}

get_as_estimator <- function(mle_estimator, c_as) {
  mle_estimator * c_as
}

get_ec_estimator <- function(mle_estimator, alpha_star) {
  mle_estimator / alpha_star
}