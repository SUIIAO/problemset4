
## functions from part a of PS3 Q2
z_test = function(X, beta, sigma = 1, mc_rep = 100){
  # function to test whether beta is equal to 0 by z-test
  #
  # Args:
  #  X: a n*p matrix has the distribution N(0, SIGMA)
  #  beta: the true coefficients of X by p and 1 matrix
  #  sigma: the error variance of Y and the default value is 1
  #  mc_rep: the number of Monte Carlo replicates and here mc_rep = 100
  #  
  # Details: 
  #  use the norm to create sigma for mc_rep times to test whether z-test
  #  precise enough to judge whether beta = 0.
  #
  # Returns: 
  #  a p-value matrix from beta_1 to beta_p with mc_rep times different result.
  
  # Get basic data of X and beta : --------------------------------------------
  n = dim(X)[1]
  p = dim(X)[2]
  
  p_value = rep(0, p * mc_rep)
  dim(p_value) = c(p, mc_rep)
  for (i in 1:mc_rep) {
    # Compute y: ----------------------------------------------------------------
    Y = X %*% beta + sigma * rnorm(n)
    
    # Use QR Decomposition to calculate the beta_hat: ---------------------------
    QR = qr(t(X) %*% X)
    R = qr.R(QR)
    Q = qr.Q(QR)
    beta_hat = solve(qr.R(QR), t(qr.Q(QR)) %*% t(X) %*% Y)
    
    # Use Y and Y_hat=X*¦Â_hat to estimate the error variance for each Monte Carlo
    # trial m
    sigma_hat = sum((Y - X %*% beta_hat)^2) / (n - p)
    
    # Calculate the p-value: ------------------------------------------------------
    v = sigma_hat * diag(chol2inv(chol(t(X) %*% X)))
    z = beta_hat / sqrt(v)
    p_value[,i] = 2 * (1 - pnorm(abs(z)))
  }
  return(p_value)
}


##  functions from part c of PS3 Q2
evaluate = function(p_value, indices = NULL){
  # function to evaluate the results of p-value and calculate their family 
  # wise error rate, false discovery rate, sensitivity and specificity
  #
  # Args:
  #  p_value: a p*mc_rep matrix carries the result of z-test
  #  indices: a set of indices where beta != 0 
  #  beta: the true coefficients of X by p and 1 matrix
  #
  # Details: 
  #  Use the Classification of multiple hypothesis tests to calculate 
  #  their family wise error rate, false dscovery rate, sensitivity and
  #  specificity.
  #
  # Returns: 
  #  a evaluation data.table carries those four value
  
  p = dim(p_value)[1]
  mc_rep = dim(p_value)[2]
  
  # Prepare a Null evaluation data.table to store the Classification of 
  # multiple hypothesis tests
  evaluation = data.table(TN = numeric(), FP = numeric(),
                          FN = numeric(), TP = numeric())
  for (i in 1:mc_rep) {
    TN = sum(p_value[-indices,i] >= 0.05)
    FP = sum(p_value[-indices,i] <  0.05)
    FN = sum(p_value[indices, i] >= 0.05)
    TP = sum(p_value[indices, i] <  0.05)
    evaluation = rbind(evaluation, list(TN, FP, FN, TP))
  }
  
  # Calculate the family wise error rate, false discovery rate, sensitivity and
  # specificity.
  evaluation = evaluation %>% mutate("Family wise error rate" = 1 * (FP >= 1),
                                     "False discovery rate" = ifelse(
                                       is.nan(FP / (FP + TP)), 0, 
                                       FP / (FP + TP)),
                                     Sensitivity = TP / (TP + FN),
                                     Specificity = TN / (TN + FP)) %>%
    summarise(Family_se = sd(`Family wise error rate`)/sqrt(mc_rep),
              False_se = sd(`False discovery rate`)/sqrt(mc_rep),
              Sensitivity_se = sd(Sensitivity)/sqrt(mc_rep),
              Specificity_se = sd(Specificity)/sqrt(mc_rep),
              "Family wise error rate" = mean(`Family wise error rate`),
              "False discovery rate" = mean(`False discovery rate`), 
              Sensitivity = mean(Sensitivity),
              Specificity = mean(Specificity))
 
  # Reorganize the results into a long data frame with columns:
  # "metric", "est", and "se"
  est = evaluation %>% select(`Family wise error rate`:Specificity) %>%
    gather(value = "est", key = "metric",
           `Family wise error rate`:Specificity)
  se = evaluation %>% select(Family_se:Specificity_se) %>%
    gather(value = "se", key = "metric_se",
           Family_se:Specificity_se)
  evaluation = cbind(est, se = se$se)
  return(evaluation)
}

method = function(p_value, indices = NULL){
  # function to evaluate the results of p-value and calculate their family 
  # wise error rate, false discovery rate, sensitivity and specificity by
  # four multiple comparison methods (bonferroni, holm, Benjamini-Hochberg
  # Benjamini-Yekuteli)
  #
  # Args:
  #  p_value: a p*mc_rep matrix carries the result of z-test
  #  indices: a set of indices where beta != 0 
  #  beta: the true coefficients of X by p and 1 matrix
  #
  # Details: 
  #  Use the command p.adjust to adjust p-values by Bonferroni, Holm, 
  #  Benjamini-Hochberg, Benjamini-Yekuteli and then use function evaluate
  #  to calculate their family wise error rate, false dscovery rate, 
  #  sensitivity and specificity.
  #
  # Returns: 
  #  a result data.table carries those values
  
  # Use bonferroni multiple comparison to correct the p-value
  p_value_bon = apply(p_value, 2, function(p) p.adjust(p, method = "bonferroni"))
  bon_e = evaluate(p_value_bon, indices) 
  bon_e = data.table(method = "Bonferroni",bon_e)
  
  # Use Holm multiple comparison to correct the p-value
  p_value_holm = apply(p_value, 2, function(p) p.adjust(p, method = "holm"))
  holm_e = evaluate(p_value_holm, indices)
  holm_e = data.table(method = "Holm",holm_e)
  
  # Use Benjamini-Hochberg multiple comparison to correct the p-value
  p_value_bh = apply(p_value, 2, function(p) p.adjust(p, method = "BH"))
  bh_e = evaluate(p_value_bh, indices)
  bh_e = data.table(method = "Benjamini-Hochberg",bh_e)
  
  # Use Benjamini-Yekuteli multiple comparison to correct the p-value
  p_value_by = apply(p_value, 2, function(p) p.adjust(p, method = "BY"))
  by_e = evaluate(p_value_by, indices)
  by_e = data.table(method = "Benjamini-Yekuteli",by_e)
  
  # Combine the result into one data.table
  result = rbind(bon_e, holm_e)
  result = rbind(result, bh_e)
  result = rbind(result, by_e)
  
  return(result)
}