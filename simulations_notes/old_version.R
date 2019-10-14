# data: 100x10
# petla: 1000

alfa_star <- 1.1678
sigma_star <- 3.3466
ile <- 100

library(sigmoid)

behaviour <- function()
{
  n <- 100 #liczba obserwacji(prob)
  kappa <- 0.1 #p/n
  p <- round(kappa*n) #liczba zmiennych
  eps <- 0.5 #czesc zmiennych niezerowych
  k <- round(eps*p) #liczba zmiennych z niezerowymi wspolczynnikami
  gamma2 <- 5 #gamma^2
  gamma <- sqrt(gamma2) #gamma
  amplitude <- sqrt(gamma2)/(sqrt(kappa*eps)) #wartosci niezerowych wspolczynnikow
  
  beta <- c(rep(amplitude,times = k), rep(0,times = p - k)) #wektor wspolczynnikow, p-k z nich to 0, a k to 'amplitude'
  data <- data.frame(matrix(0,p,ile))
  
  c_OS <- (n*gamma2*alfa_star) / ((sigma_star^2)*p + n*gamma2*(alfa_star^2))
  
  z <- list()
  
  
  for(i in 1:ile)
  {
    X <- matrix(rnorm(n*p)/sqrt(n), n, p)
    P <- sigmoid(X%*%beta)
    Y <- rbinom(n, 1, P)
    model <- glm(Y~X,family = binomial)
    beta_MLE <- model$coef[2:(p+1)] #estymowane wspolczynniki
    data[,i] <- beta_MLE 
  }
  
  biasMLE<-c() 
  mseMLE <- c() 
  newMLE <- c() 
  biasEC<-c() 
  mseEC <- c()
  newEC <- c() 
  biasOS<-c() 
  mseOS <- c() 
  newOS <- c() 
  biasOS2<-c() 
  mseOS2 <- c() 
  newOS2 <- c() 
  
  for(i in 1:ile)
  {
    biasMLE[i]<- sum(data[,i] - beta)/p
    mseMLE[i]<-sum((data[,i] - beta)^2)/p
    newMLE[i]<- sum((data[,i] - beta)*beta)/p
    
    biasEC[i]<- sum(data[,i]/alfa_star - beta)/p
    mseEC[i]<-sum((data[,i]/alfa_star - beta)^2)/p
    newEC[i]<- sum((data[,i]/alfa_star - beta)*beta)/p
    
    biasOS[i]<- sum(data[,i]*c_OS - beta)/p
    mseOS[i]<-sum((data[,i]*c_OS - beta)^2)/p
    newOS[i]<- sum((data[,i]*c_OS - beta)*beta)/p
    
  }
  
  z$biasEC <- biasEC
  z$mseEC <- mseEC
  z$newEC <-  newEC
  z$biasMLE <- biasMLE
  z$mseMLE <- mseMLE
  z$newMLE <-  newMLE
  z$biasOS <- biasOS
  z$mseOS <- mseOS
  z$newOS <-  newOS
  z$data <- data
  
  return(z)
}
