
## functions from part a of PS3 Q2
z_test = function(X, beta, sigma = 1, mc_rep = 100){
  # Simulate Y from Y|X ~ N(XB, sigma^2 I) and compute p-values corresponding to
  # Wald tests for B != 0. Repeat mc_rep times.
  #
  # Arguments:
  #   X : an n by p numeric matrix
  #   beta: a p by 1 numeric matrix
  #   sigma: std deviation for Y|X,  Y|X ~ N(XB, sigma^2 I)
  #   mc_rep: The number of Monte Carlo replications to use
  #
  # Output: A p by mc_rep matrix of p-values
  
  # This part doesn't need to change for each replication
  QR = qr( crossprod(X) )
  QX = X %*% qr.Q(QR) 
  XtXinv = solve( qr.R(QR), t( qr.Q(QR) ))
  
  n = nrow(X)
  p = ncol(X)
  
  # Generate mc_rep copies of Y at once, each in a column.
  Y = as.numeric(X %*% beta) + rnorm(n*mc_rep)
  dim(Y) = c(n, mc_rep)
  
  # estimate betas and residual standard errors
  b = solve(qr.R(QR), crossprod( QX, Y ) )
  
  # It's okay if you divide by {n - p} outside the sum, but this
  # is more comparable to what is done by .lm.fit()
  s_sq = colSums( {Y - as.numeric(X %*% b)}^2 / {n - p})
  
  # standard error of b
  v = sqrt( diag(XtXinv) * rep(s_sq, each = p) )
  
  # return a matirx of p-values
  # Use pt to replicate lm, but the normal approximation is fine here. 
  matrix( 2*pt( abs( b / v ), df = {n-p}, lower.tail = FALSE ), p, mc_rep )  
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