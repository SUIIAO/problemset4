## Problem Set 4, Question 2a
##
## Author: SU I IAO
## Updated: December 3rd, 2018

# libraries: ------------------------------------------------------------------
library(data.table)
library(dplyr)
library(tidyverse)
library(parallel)


# Question a: Use mclapply to run parallel simulations for 
# rho %in% {0.25*{-3:3}}. Let sigma = 1 and use 10,000 Monte Carlo replications
# Reorganize the results into a long data frame results_q4a with columns: 
# "rho", "sigma", "metric", "method", "est", and "se". 

# sources ps4_q2_funcs.R to get function from part a and c of PS3 Q2
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
p_value = mclapply(parameter, z_test, beta = beta, sigma = 1, mc_rep = 10000)

# Evaluate the results of p-value and calculate their family wise error rate,
# false discovery rate, sensitivity and specificity by four multiple comparison
# methods (bonferroni, holm, Benjamini-Hochberg, Benjamini-Yekuteli) 
evaluation = mclapply(p_value, method, indices = which(beta != 0))

# Reorganize the results into a long data frame results_q4a with columns:
# "rho", "sigma", "metric", "method", "est", and "se".
results_q4a = data.table()
for (i in c(-3:3)) {
  a = data.table(rho = i*0.25, sigma = 1, evaluation[[i+4]])
  if(length(results_q4a) == 0)
    results_q4a = a
  else
    results_q4a = rbind(results_q4a, a)
}
