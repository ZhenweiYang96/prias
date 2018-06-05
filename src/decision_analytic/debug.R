
hazardFunc = function (s, patientId) {
  weibullShape = 6
  weibullScale = 4
  
  b_subject = b[patientId, ]
  Age = simDs.id$Age[patientId]  
  wGamma = simDs.id$wGamma[patientId]
  
  # pdf_s = dweibull(s, shape = weibullShape, scale = weibullScale)
  # survival_s = (1-pweibull(q = s, shape= weibullShape, scale = weibullScale))
  # baselinehazard_s = pdf_s/survival_s
  
  baselinehazard_s = (weibullShape/weibullScale)*(s/weibullScale)^(weibullShape-1)
  #baselinehazard_s = exp(splineDesign(mvJoint_dre_psa_dre_value$control$knots, s, 
  #                 ord = mvJoint_dre_psa_dre_value$control$ordSpline, 
  #                 outer.ok = T) %*% mvJoint_dre_psa_dre_value$statistics$postMeans$Bs_gammas)
  
  df_s = data.frame(visitTimeYears = s, Age = Age)
  
  xi_s_logodds_high_dre_val = model.matrix(fixedDREFormula, df_s)
  xi_s_log2psaplus1_val = model.matrix(fixedPSAFormula, df_s)
  xi_s_log2psaplus1_slope = model.matrix(fixedPSASlopeFormula, df_s)
  
  zi_s_logodds_high_dre_val = model.matrix(randomDREFormula, df_s)
  zi_s_log2psaplus1_val = model.matrix(randomPSAFormula, df_s)
  zi_s_log2psaplus1_slope = model.matrix(randomPSASlopeFormula, df_s)
  
  #There are 7 random effects, 1,2 for DRE and 3,4,5,6,7 for PSA
  zib_logodds_high_dre_val = zi_s_logodds_high_dre_val %*% b_subject[1:2]
  zib_log2psaplus1_val = zi_s_log2psaplus1_val %*% b_subject[3:7]
  zib_log2psaplus1_slope = zi_s_log2psaplus1_slope %*% b_subject[4:7] #One less random effect to ignore intercept
  
  xBetaZb_s_logodds_high_dre_value = xi_s_logodds_high_dre_val  %*% betas_dre + zib_logodds_high_dre_val
  xBetaZb_s_log2psaplus1_value = xi_s_log2psaplus1_val %*% betas_psa + zib_log2psaplus1_val
  xBetaZb_s_log2psaplus1_slope = xi_s_log2psaplus1_slope %*% betas_psa[-c(1:3)] + zib_log2psaplus1_slope #-c(1:3) to ignore intercept, age and age^2
  
  y_Alpha = cbind(xBetaZb_s_logodds_high_dre_value, xBetaZb_s_log2psaplus1_value, xBetaZb_s_log2psaplus1_slope) %*% alphas
  #y_Alpha = cbind(xBetaZb_s_value) %*% alphas[1]
  
  baselinehazard_s * exp(wGamma + y_Alpha)
}