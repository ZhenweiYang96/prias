load("Rdata/gap3/PRIAS_2019/cleandata.Rdata")

fixed_random_psaSlopeFormula = ~ 0 + dns(I(year_visit-2)/2, knots=(c(0.5, 1.3, 3)-2)/2, Boundary.knots=(c(0, 6.3)-2)/2)

prias_psa = prias_long_final[!is.na(prias_long_final$psa),]
X_Z = model.matrix(fixed_random_psaSlopeFormula, data = prias_psa)
psacount_per_patient = table(prias_psa$P_ID)
fitted_velocities = X_Z %*% mvJoint_psa_time_scaled$statistics$postMeans$betas1[-c(1:2)] +
  apply(X_Z * mvJoint_psa_time_scaled$statistics$postMeans$b[rep(1:7813, psacount_per_patient),-1], MARGIN = 1, sum)
summary(fitted_velocities)  
#round(quantile(fitted_velocities, probs = c(0.025,0.975)),3)
