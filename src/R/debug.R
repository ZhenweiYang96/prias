hazardFunc = function (s, i) {
  pdf_s = dweibull(s, shape = weibullShape, scale = weibullScale)
  survival_s = (1-pweibull(q = s,shape= weibullShape, scale = weibullScale))
  baselinehazard_s = pdf_s/survival_s

  df_s = data.frame(visitTimeYears = s, Age = simDs.id$Age[i])
  
  xi_s_val = model.matrix(fixedValueFormula, df_s)
  xi_s_slope = model.matrix(fixedSlopeFormula, df_s)
  
  zi_s_val = model.matrix(randomValueFormula, df_s)
  zi_s_slope = model.matrix(randomSlopeFormula, df_s)
  
  zib_val = zi_s_val %*% b[i, ]
  zib_slope = zi_s_slope %*% b[i, -1] #-1 to ignore intercept
  
  xBetaZb_s_value = xi_s_val %*% betas + zib_val
  xBetaZb_s_slope = xi_s_slope %*% betas[-c(1:3)] + zib_slope #-c(1:3) to ignore intercept, age and age^2
  
  y_Alpha = cbind(xBetaZb_s_value, xBetaZb_s_slope) %*% getAlpha(fittedJointModel, weightedOnes = F)
  #y_Alpha = cbind(xBetaZb_s_value) %*% getAlpha(fittedJointModel, weightedOnes = F)
  
  c(baselinehazard_s * exp(wGamma[i] + y_Alpha))
}