## Problem Set 4, Question 2b
##
## Author: SU I IAO
## Updated: December 3rd, 2018

# libraries: ------------------------------------------------------------------
library(data.table)
library(dplyr)
library(tidyverse)
library(parallel)
library(doParallel)

# Question b: Setup a 4 core cluster using doParallel and then use nested 
# foreach loops to run simulations for rho %in% {0.25*{-3:3}} and 
# sigma = {.25,.5,1}. Reshape the results as before into results_q4b 
# saved to a file  results_q4b.RData. 
# Use a PBS file to run this script on the Flux cluster.

# sources ps4_q2_funcs.R to get functions from part a and c of PS3 Q2 and 
# function "method"
source('ps4_q2_funcs.R')

# Setup a 4 core cluster using doParallel
ncores = 4  
cl = makeCluster(ncores)
registerDoParallel(cl)

## do parallel computaitons with foreach
results = foreach(sigma = c(0.25, 0.5, 1)) %dopar% {
  library(data.table)
  library(dplyr)
  library(tidyverse)
  library(parallel)
  
  source('ps4_q2_funcs.R')
  # Build x_sigma by given rho and then generate X randomly by x_sigma: ---------
  n = 1000
  p = 100
  beta = c(rep(.1, 10), rep(0, p - 10))
  parameter = list()
  x_sigma = rep(0, 100*100)
  dim(x_sigma) = c(100, 100)
  diag(x_sigma) = rep(1, 100)
  
  for (i in c(-3:3)) {
    rho = 0.25*i
    SIGMA_10 = matrix(rep(0.01*rho,100), 10, 10)
    diag(SIGMA_10) = rep(1, 10)
    x_sigma[1:10, 1:10] = SIGMA_10
    R = chol(x_sigma)
    X = rnorm(n * p)
    dim(X) = c(n, p)
    X = X %*% R
    parameter = c(parameter, list(X))
  }
  
  # Use mclapply and function z_test to get p-value under different rho
  p_value = mclapply(parameter, z_test, beta = beta, sigma = sigma, mc_rep = 10000)
  
  # Evaluate the results of p-value and calculate their family wise error rate,
  # false discovery rate, sensitivity and specificity by four multiple comparison
  # methods (bonferroni, holm, Benjamini-Hochberg, Benjamini-Yekuteli) 
  evaluation = mclapply(p_value, method, indices = which(beta != 0))
  
  # Reorganize the results into a long data frame results_q4a with columns:
  # "rho", "sigma", "metric", "method", "est", and "se".
  results_q4a = data.table()
  for (i in c(-3:3)) {
    a = data.table(rho = i*0.25, sigma = sigma, evaluation[[i+4]])
    if(length(results_q4a) == 0)
      results_q4a = a
    else
      results_q4a = rbind(results_q4a, a)
  }
  results_q4a
}
stopCluster(cl)

results_q4b = rbind(results[[1]], results[[2]])
results_q4b = rbind(results_q4b, results[[3]])
save(results_q4b, file="results_q4b.Rdata")
     