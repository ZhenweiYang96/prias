library(JMbayes)
load("Rdata/gap3/PRIAS_2019/cleandata.Rdata")
load("Rdata/lastpaper/fitted_model/mvJoint_dre_psa_2knots_quad_age_light.Rdata")

fixed_random_psaSlopeFormula = ~ 0 + dns(year_visit, knots=c(0.75, 2.12), Boundary.knots=c(0, 6.4))

prias_psa = prias_long_final[!is.na(prias_long_final$psa),]
X_Z = model.matrix(fixed_random_psaSlopeFormula, data = prias_psa)
psacount_per_patient = table(prias_psa$P_ID)
fitted_velocities = X_Z %*% mvJoint_dre_psa_2knots_lin_age$statistics$postMeans$betas2[-c(1:2)] +
  apply(X_Z * mvJoint_dre_psa_2knots_lin_age$statistics$postMeans$b[rep(1:7813, psacount_per_patient),-c(1:3)], MARGIN = 1, sum)
summary(fitted_velocities)  
#round(quantile(fitted_velocities, probs = c(0.025,0.975)),3)
