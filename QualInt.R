library(mvtnorm)
library(expm)

NaiveTest <- function(theta.hat, standard.error, alpha) {
  sep.p.vals <- pchisq(theta.hat^2/standard.error^2, df = 1, lower.tail = FALSE)
  p.val <- ifelse(((sep.p.vals[1] <  alpha) & (sep.p.vals[2] > alpha)) |
                  ((sep.p.vals[1] >  alpha) & (sep.p.vals[2] < alpha)),
                    0, 1)
  # Test does not provide real p-values
  # We set p.val = 0 for "reject" and p.val = 1 for "fail to reject"
  return(p.val)
}

RelDifLRTStat <- function(theta.hat, standard.error, kappa = 2) {
  
  theta.hat.1 <- abs(theta.hat[1])
  theta.hat.2 <- abs(theta.hat[2])
  
  theta.hat.max <- (theta.hat.1 + theta.hat.2)/2 + abs(theta.hat.1 - theta.hat.2)/2
  theta.hat.min <- (theta.hat.1 + theta.hat.2)/2 - abs(theta.hat.1 - theta.hat.2)/2
  
  tau.hat.max <- standard.error[1] ** 2 * (theta.hat.1 > theta.hat.2) + 
                 standard.error[2] ** 2 * (theta.hat.1 < theta.hat.2)
  tau.hat.min <- standard.error[1] ** 2 * (theta.hat.1 < theta.hat.2) + 
                 standard.error[2] ** 2 * (theta.hat.1 > theta.hat.2)
  
  test.stat <- (theta.hat.max - kappa * theta.hat.min) * 
               (tau.hat.max + kappa ** 2 * tau.hat.min)^(-.5)
  
  return(test.stat)
}

RelDifLRTLocal <- function(test.stat, asymp.mean, asymp.var, lambda, kappa = 2) {
  
  c.1 <- asymp.mean[1]
  c.2 <- asymp.mean[2]
  v.1 <- asymp.var[1]
  v.2 <- asymp.var[2]
  
  cstar.11 <- ((1 - lambda) * c.1 - kappa * lambda * c.2)/
              sqrt((1 - lambda) * v.1 + kappa ** 2 * lambda * v.2)
  cstar.12 <- ((1 - lambda) * c.1 + kappa * lambda * c.2)/
              sqrt((1 - lambda) * v.1 + kappa ** 2 * lambda * v.2)
  cstar.21 <- (-kappa * (1 - lambda) * c.1 + lambda * c.2)/
              sqrt(kappa ** 2 * (1 - lambda) * v.1 + lambda * v.2)
  cstar.22 <- (kappa * (1 - lambda) * c.1 + lambda * c.2)/
              sqrt(kappa ** 2 * (1 - lambda) * v.1 + lambda * v.2)
  
  nu.1 <- ((1 - lambda) * v.1 - kappa ** 2 * lambda * v.2)/
          ((1 - lambda) * v.1 + kappa ** 2 * lambda * v.2)
  nu.2 <- (lambda * v.2 - kappa ** 2 * (1 - lambda) * v.1)/
          (kappa ** 2 * (1 - lambda) * v.1 + lambda * v.2)
  
  out <-  pmvnorm(lower = c(test.stat, test.stat), upper = c(Inf, Inf),
                  mean = c(cstar.11, cstar.12), 
                  sigma = matrix(c(1, nu.1, nu.1, 1), 2, 2)) +
          pmvnorm(lower = c(-Inf, -Inf), upper = c(-test.stat, -test.stat),
                  mean = c(cstar.11, cstar.12), 
                  sigma = matrix(c(1, nu.1, nu.1, 1), 2, 2)) +
          pmvnorm(lower = c(test.stat, test.stat), upper = c(Inf, Inf),
                  mean = c(cstar.21, cstar.22), 
                  sigma = matrix(c(1, nu.2, nu.2, 1), 2, 2)) +
          pmvnorm(lower = c(-Inf, -Inf), upper = c(-test.stat, -test.stat),
                  mean = c(cstar.21, cstar.22), 
                  sigma = matrix(c(1, nu.2, nu.2, 1), 2, 2))
  
  return(out ** (test.stat > 0))
}

RelDifLRT <- function(theta.hat, standard.error, n, kappa = 2) {
  
  test.stat <- RelDifLRTStat(theta.hat, standard.error, kappa)
  p.val <- max(c(pchisq(test.stat ** 2, df = 1, lower.tail = FALSE),
                 RelDifLRTLocal(test.stat, c(0,0), n * standard.error ** 2,
                                lambda = n[1]/sum(n))))
  
  out <- list(p.val = p.val,
              test.stat = test.stat)
  return(out)
}

kappa.max <- function(theta.hat, standard.error, n, alpha, tol = .Machine$double.eps^0.25) {
  
  kappa.lower <- ifelse(RelDifLRT(theta.hat, standard.error, n, 1)$p.val < alpha, 1, NA)
  kappa.upper <- 4
  while(RelDifLRT(theta.hat, standard.error, n, kappa.upper)$p.val < alpha)
    kappa.upper <- kappa.upper * 1.1
  
  if(!is.na(kappa.lower)) {
    kappa.max <- uniroot(function(x) RelDifLRT(theta.hat, standard.error, n, x)$p.val - alpha,
                         interval = c(kappa.lower, kappa.upper), tol = tol)$root
  }
  else {
    kappa.max <- 1
  }
  
  return(kappa.max)
}

OmnibusLRTStat <- function(theta.hat, standard.error, kappa = 2) {
  
  test.stat <- min(c((theta.hat[1] - kappa * theta.hat[2]) ** 2/
                     (standard.error[1] ** 2 + kappa ** 2 * standard.error[2] ** 2),
                     (kappa * theta.hat[1] - theta.hat[2]) ** 2/
                     (kappa ** 2 * standard.error[1] ** 2 + standard.error[2] ** 2)))
  test.stat <- ifelse(theta.hat[1] * theta.hat[2] < 0 | 
                      abs(theta.hat[1]/theta.hat[2]) > kappa |
                      abs(theta.hat[2]/theta.hat[1]) > kappa,
                      test.stat, 0)
  
  return(test.stat)
}

OmnibusLRTLocal <- function(test.stat, asymp.mean, asymp.var, lambda, kappa = 2) {
  
  c.1 <- asymp.mean[1]
  c.2 <- asymp.mean[2]
  v.1 <- asymp.var[1]
  v.2 <- asymp.var[2]
  
  cstar.1 <- ((1 - lambda) * c.1 - kappa * lambda * c.2)/
             sqrt((1 - lambda) * v.1 + kappa ** 2 * lambda * v.2)
  cstar.2 <- (kappa * (1 - lambda) * c.1 - lambda * c.2)/
             sqrt(kappa ** 2 * (1 - lambda) * v.1 + lambda * v.2)
  nu <- (kappa * (1 - lambda) * v.1 + kappa * lambda * v.2)/
        (sqrt((1- lambda) * v.1 + kappa ** 2 * lambda * v.2) *
         sqrt(kappa ** 2 * (1 - lambda) * v.1 + lambda * v.2))
  
  out <-  pmvnorm(lower = sqrt(c(test.stat, test.stat)), upper = c(Inf, Inf),
            mean = c(cstar.1, cstar.2), 
            sigma = matrix(c(1, nu, nu, 1), 2, 2)) +
          pmvnorm(lower = c(-Inf, -Inf), upper = -sqrt(c(test.stat, test.stat)),
                  mean = c(cstar.1, cstar.2), 
                  sigma = matrix(c(1, nu, nu, 1), 2, 2))
  
  return(out)
}

OmnibusLRT <- function(theta.hat, standard.error, n, kappa = 2) {
  
  test.stat <- OmnibusLRTStat(theta.hat, standard.error, kappa)
  p.val <- max(c(pchisq(test.stat ** 2, df = 1, lower.tail = FALSE),
                 OmnibusLRTLocal(test.stat, c(0,0), n * standard.error ** 2,
                                 n[1]/sum(n))))
  
  out <- list(p.val = p.val,
              test.stat = test.stat)
  return(out)
}
