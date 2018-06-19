makeTimeDepDs = function(simDs.id){
  simDs.id = simDs.id[order(simDs.id$progression_time, decreasing = F),]
  
  simDs_time_dep = simDs.id[rep(1:nrow(simDs.id), 1:nrow(simDs.id)),]
  simDs_time_dep$startTime = c(0,simDs.id$progression_time[-nrow(simDs.id)])[unlist(sapply(1:nrow(simDs.id), function(x){1:x}))]
  simDs_time_dep$visitTimeYears = simDs.id$progression_time[unlist(sapply(1:nrow(simDs.id), function(x){1:x}))]
  simDs_time_dep$progressed_in_period = unlist(sapply(1:nrow(simDs.id), function(x){c(rep(0, x-1),simDs.id$progressed[x])}))
  
  return(simDs_time_dep)
}    

#When you call this function replace the contents of simDs_time_dep colnames c("b_Int_PSA", "b_Slope1_PSA", "b_Slope2_PSA")
getTimeDepDsCoefficients = function(simDs_time_dep, betas){
  fixedPSAFormula = ~ 1 + I(Age - 70) +  I((Age - 70)^2) + ns(visitTimeYears, knots=c(0.1, 0.5, 4), Boundary.knots=c(0, 7))
  randomPSAFormula = ~ 1 + ns(visitTimeYears, knots=c(0.1), Boundary.knots=c(0, 7))
  
  fixedPSASlopeFormula = ~ 0 + dns(visitTimeYears, knots=c(0.1, 0.5, 4), Boundary.knots=c(0, 7))
  randomPSASlopeFormula = ~ 0 + dns(visitTimeYears, knots=c(0.1), Boundary.knots=c(0, 7))
  
  simDs_time_dep$true_PSA_in_period = model.matrix(fixedPSAFormula, simDs_time_dep) %*% betas + rowSums(model.matrix(randomPSAFormula, simDs_time_dep) * simDs_time_dep[, c("b_Int_PSA", "b_Slope1_PSA", "b_Slope2_PSA")])
  simDs_time_dep$true_PSASlope_in_period = model.matrix(fixedPSASlopeFormula, simDs_time_dep) %*% betas[4:7] + rowSums(model.matrix(randomPSASlopeFormula, simDs_time_dep) * simDs_time_dep[, c("b_Slope1_PSA", "b_Slope2_PSA")])
  
  relative_risk_fitted  = coxph(Surv(startTime, visitTimeYears, progressed_in_period) ~ I(Age - 70) +  I((Age - 70)^2) + 
                                  true_PSA_in_period + true_PSASlope_in_period, data=simDs_time_dep)
  return(relative_risk_fitted$coefficients)
}
