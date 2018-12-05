## Problem Set 4, Question 2a
##
## Author: SU I IAO
## Updated: December 3rd, 2018

# libraries: ------------------------------------------------------------------
library(data.table)
library(dplyr)
library(tidyverse)
library(parallel)
library(future)

# Question c: Modify the script from part a to create ps4_q2c.R which reads the
# following arguments from the command line: sigma, mc_rep, and  n_cores. 
# Also modify the script to use the futures package for parallelism. 
# Use a PBS file to run this script as a job array for sigma={.25,.5,1}.

# default arguments
args_list = list(
  sigma = 1,
  mc_rep = 10000,
  n_cores = 4)

## get parameters from command line
args = commandArgs(trailingOnly = TRUE)

# functions for finding named arguments
args_to_list = function(args){
  ind = grep('=', args)  
  args_list = sapply(args[ind],strsplit, '=')
  names(args_list) = sapply(args_list, function(x) x[1])
  
  args_list = lapply(args_list, function(x) as.numeric(x[2]))
  args_list
}

# get named arguments
args_list_in = args_to_list(args)

# update non default arguments
ignored = c()
for ( arg in names(args_list_in) ) {
  # Check for unknown argument
  if ( is.null(args_list[[arg]]) ) {
    ignored = c(ignored, arg)
  } else{
    # update if known
    args_list[[arg]] = args_list_in[[arg]]
  }
}
sigma = args_list[[1]]
mc_rep = args_list[[2]]
n_cores = args_list[[3]]

# sources ps4_q2_funcs.R to get function from part a and c of PS3 Q2
source('ps4_q2_funcs.R')

# Build x_sigma by given rho and then generate X randomly by x_sigma: ---------
n = 1000
p = 100
beta = c(rep(.1, 10), rep(0, p - 10))
x_sigma = rep(0, 100*100)
dim(x_sigma) = c(100, 100)
diag(x_sigma) = rep(1, 100)

plan(multisession)
parameter = list()
for (i in c(-3:3)) {
  parameter[[i+4]] = future({
    rho = 0.25*i
    SIGMA_10 = matrix(rep(0.01*rho,100), 10, 10)
    diag(SIGMA_10) = rep(1, 10)
    x_sigma[1:10, 1:10] = SIGMA_10
    R = chol(x_sigma)
    X = rnorm(n * p)
    dim(X) = c(n, p)
    X = X %*% R
    X
  })
}
parameter = lapply(parameter, value)
# Use mclapply and function z_test to get p-value under different rho
p_value = mclapply(parameter, z_test, beta = beta, sigma = sigma, mc_rep = mc_rep)

# Evaluate the results of p-value and calculate their family wise error rate,
# false discovery rate, sensitivity and specificity by four multiple comparison
# methods (bonferroni, holm, Benjamini-Hochberg, Benjamini-Yekuteli) 
evaluation = mclapply(p_value, method, indices = which(beta != 0))

# Reorganize the results into a long data frame results_q4a with columns:
# "rho", "sigma", "metric", "method", "est", and "se".
results_q4c = data.table()
for (i in c(-3:3)) {
  c = data.table(rho = i*0.25, sigma = sigma, evaluation[[i+4]])
  if(length(results_q4a) == 0)
    results_q4c = c
  else
    results_q4c = rbind(results_q4c, c)
}
