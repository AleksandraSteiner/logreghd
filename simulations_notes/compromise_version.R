# data: 100x10
# petla: 1000

library(sigmoid)

posrednia <- function(ile = 100, n = 100, kappa = 0.1, eps = 0.5) {
  alfa_star <- 1.1678
  sigma_star <- 3.3466
  
  p <- round(kappa*n) #liczba zmiennych
  k <- round(eps*p) #liczba zmiennych z niezerowymi wspolczynnikami
  gamma2 <- 5 #gamma^2
  gamma <- sqrt(gamma2) #gamma
  amplitude <- sqrt(gamma2)/(sqrt(kappa*eps)) #wartosci niezerowych wspolczynnikow
  
  beta <- c(rep(amplitude,times = k), rep(0,times = p - k)) #wektor wspolczynnikow, p-k z nich to 0, a k to 'amplitude'
  data <- data.frame(matrix(0,p,ile))
  
  c_OS <- (n*gamma2*alfa_star) / ((sigma_star^2)*p + n*gamma2*(alfa_star^2))
  
  data <- apply(data, 2, function(i) {
    X <- matrix(rnorm(n*p)/sqrt(n), n, p)
    P <- sigmoid(X%*%beta)
    Y <- rbinom(n, 1, P)
    model <- glm(Y~X,family = binomial)
    model$coef[2:(p+1)] #estymowane wspolczynniki
  })
  
  lapply(data, function(estimated_beta) {
    list(
      biasMLE = sum(estimated_beta - beta)/p,
      mseMLE = sum((estimated_beta - beta)^2)/p,
      newMLE = sum((estimated_beta - beta)*beta)/p,
      biasEC = sum(estimated_beta/alfa_star - beta)/p,
      mseEC = sum((estimated_beta/alfa_star - beta)^2)/p,
      newEC = sum((estimated_beta/alfa_star - beta)*beta)/p,
      biasOS = sum(estimated_beta*c_OS - beta)/p,
      mseOS = sum((estimated_beta*c_OS - beta)^2)/p,
      newOS = sum((estimated_beta*c_OS - beta)*beta)/p
    )
  })
}



