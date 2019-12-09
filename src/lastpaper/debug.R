seed = c(2021:2030, 2041:2044)

year_visit = seq(0.1, 10, length.out = 100)

load("Rdata/lastpaper/fitted_model/mvJoint_dre_psa_2knots_quad_age_light.Rdata")

baselinehazard = exp(splineDesign(mvJoint_dre_psa_2knots_quad_age_light$control$knots, year_visit, 
                                  ord = mvJoint_dre_psa_2knots_quad_age_light$control$ordSpline, outer.ok = T) %*% mvJoint_dre_psa_2knots_quad_age_light$statistics$postMeans$Bs_gammas)

tt = sapply(seed, FUN = function(seed){
  load(paste0("Rdata/lastpaper/simulation/light/jointModelData_seed_", seed, "_t3.Rdata"))
  baselinehazard = exp(splineDesign(jointModelData$mvJoint_dre_psa_simDs$control$knots, year_visit, 
                                    ord = jointModelData$mvJoint_dre_psa_simDs$control$ordSpline, outer.ok = T) %*% jointModelData$mvJoint_dre_psa_simDs$statistics$postMeans$Bs_gammas)
  return(baselinehazard)
})

library(ggplot2)

ggplot() + geom_line(aes(x=year_visit, y=baselinehazard), color='red') +
  geom_line(aes(x=year_visit, y=rowMeans(tt))) + 
  geom_ribbon(aes(x=year_visit, ymin=apply(tt,1, quantile, probs=0.025),
                  ymax=apply(tt,1,quantile, probs=0.975)), alpha=0.15)