set.seed(206)

my.wd <- "/Users/awhudson/Dropbox/MethodsResearch/QualitativeInteractions/Code/"
setwd(my.wd)
source('Method/RelativeDifferenceTest-rev-su22.R')

##########################################
### Simulation Settings
##########################################

theta.1 <- 1
theta.2 <- seq(-1, 1, .05)
n <- c(50, 100, 200, 400, 800, 1600)
kappa <- c(1.5, 2, 4, 8)
sigma <- 1
no.reps <- 1000
alpha <- .05

###########################################
### Obtain simulation output
###########################################

### For storing results
p.vals.RD <- array(NA, dim = c(length(theta.2), length(n), no.reps, length(kappa)))
p.vals.OM <- array(NA, dim = c(length(theta.2), length(n), no.reps, length(kappa)))
p.vals.NAIVE <- array(NA, dim = c(length(theta.2), length(n), no.reps))
p.vals.QUANT <- array(NA, dim = c(length(theta.2), length(n), no.reps))
kmax <- array(NA, dim = c(length(theta.2), length(n), no.reps))

for(i in 1:length(theta.2)) {
  print(paste0("theta.2 = ", theta.2[i]))
  for(j in 1:length(n)) {
    for(k in 1:no.reps){
      # Generate datasets for each population
      x.1 <- rnorm(n[j])
      x.2 <- rnorm(n[j])
      y.1 <- theta.1 * x.1 + rnorm(n[j], sd = sigma)
      y.2 <- theta.2[i] * x.2 + rnorm(n[j], sd = sigma)
      
      # Obtain effect estimates and standard error
      mod.1 <- lm(y.1 ~ x.1)
      mod.2 <- lm(y.2 ~ x.2)
      theta.hat <- c(coefficients(mod.1)[2],
                     coefficients(mod.2)[2])
      standard.error <- c(vcov(mod.1)[2,2] ** .5,
                          vcov(mod.2)[2,2] ** .5)
      
      # Perform test for quantitative interactions
      quant.test.stat <- (theta.hat[1] - theta.hat[2])^2/(sum(standard.error^2))
      p.vals.QUANT[i,j,k] <- pchisq(quant.test.stat, df = 1, lower.tail = FALSE)
      
      # Peform naive test
      p.vals.NAIVE[i,j,k] <- NaiveTest(theta.hat, standard.error, .05)
      
      # Perform likelihood ratio test
      for(l in 1:length(kappa)) {
        p.vals.RD[i,j,k,l] <- RelDifLRT(theta.hat, standard.error,
                                        c(n[j], n[j]), kappa[l])$p.val
        p.vals.OM[i,j,k,l] <- OmnibusLRT(theta.hat, standard.error,
                                         c(n[j], n[j]), kappa[l])$p.val
      }
      # Obtain kappa_max statistic
      kmax[i,j,k] <- kappa.max(theta.hat, standard.error,
                               c(n[j], n[j]), alpha, tol = .01)
    }
  }
}

output = list(p.vals.RD = p.vals.RD, p.vals.OM = p.vals.OM, p.vals.NAIVE = p.vals.NAIVE, kmax = kmax,
              n = n, kappa = kappa, theta.2 = theta.2)

#######################################
### Save Output
#######################################
output = list(p.vals.RD = p.vals.RD, p.vals.OM = p.vals.OM, 
              p.vals.NAIVE = p.vals.NAIVE, p.vals.QUANT = p.vals.QUANT,
              kmax = kmax, n = n, kappa = kappa, theta.2 = theta.2)
save(output,
     file = "/Users/awhudson/Dropbox/MethodsResearch/QualitativeInteractions/Code/Sims-rev-su2022/sim_output.rda")
