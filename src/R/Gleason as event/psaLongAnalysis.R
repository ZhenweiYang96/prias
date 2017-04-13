source("src/R/common.R")

psa_data_set =  prias_long[!is.na(prias_long$psa),]
psa_data_set$logpsa1 = log(psa_data_set$psa + 1, base = 2)
save.image(file = "Rdata/Gleason as event/psa_gl.Rdata")

psa_data_set$residuals = rep(NA, nrow(psa_data_set))
psa_data_set$residuals[!is.na(psa_data_set$psa)] = residuals(lm(log(psa + 1)~poly(Age,2), data = psa_data_set))

pidList = unique(psa_data_set$P_ID)
randomSample = sample(x = pidList, replace = F, size = 150)
ggplot(data=psa_data_set[psa_data_set$P_ID %in% randomSample,], aes(x=visitTimeYears, y=log(psa+1, base = 2))) + 
  geom_line(aes(group=P_ID)) + stat_smooth() + 
  ticksX(from=0,max(psa_data_set$visitTimeYears), 0.25) + xlab("Time (years)") + 
  ylab("log(psa + 1): Adjusted for Age")

#Since every box plot has 1 entry per person, 
#if the first few box plots have too many time points it means that it is not because
#people came frequently, it could also be due to the fact that people started dropping out later
#and so initially we see more concentration of time points.
#Per person the concentration is not more in the beginning, but overall marginally
#there is more concentration in the beginning
ggplot(data=psa_data_set[psa_data_set$P_ID,], aes(factor(visit_number),visitTimeYears)) + 
  geom_boxplot() + stat_summary(fun.data = function(x){
    return(c(y = 1.2, label = length(x))) 
  }, geom = "text", fun.y = median) + ticksY(max(psa_data_set$visitTimeYears), 0.05)

quantiles_time_years = quantile(psa_data_set$visitTimeYears, probs = c(0.25, 0.5, 0.75))
# 25%        50%        75% 
# 0.04849315 0.12136986 0.24849315 

###########################
# Linear evolutions
###########################
psaModel_linear = lme(fixed=logpsa1~I(Age - 70) +  I((Age - 70)^2) + 
                        visitTimeYears,
                     random = ~visitTimeYears|P_ID, 
                     data=psa_data_set,
                     control = lmeControl(opt = "optim", optimMethod = "L-BFGS-B"))

joint_psa_linear = jointModelBayes(psaModel_linear, coxModel, timeVar = "visitTimeYears",
                            param = "td-both",
                            extraForm = list(fixed = ~ 1, random= ~ 1,
                                            indFixed = 4, indRandom=2))


mvglm_psa_linear = mvglmer(list(
  logpsa1 ~  I(Age - 70) +  I((Age - 70)^2) + visitTimeYears + 
    (visitTimeYears|P_ID)), 
  data =psa_data_set, families = list(gaussian))
save.image(file = "Rdata/Gleason as event/psa_gl_linear.Rdata")

mvJoint_psa_linear_tdval = mvJointModelBayes(mvglm_psa_linear, coxModel, timeVar = "visitTimeYears")
save.image(file = "Rdata/Gleason as event/psa_gl_linear.Rdata")

mvJoint_psa_linear_tdboth = mvJointModelBayes(mvglm_psa_linear, coxModel, timeVar = "visitTimeYears", 
                                             Formulas = list("logpsa1" = "value",
                                                             "logpsa1" = list(fixed = ~ 1,
                                                                              random= ~ 1,
                                                                              indFixed = 4, indRandom=2, 
                                                                              name = "slope")))

save.image(file = "Rdata/Gleason as event/psa_gl_linear.Rdata")
################################
# Non linear evolutions
################################
psaModel_spline_default = fitUnivaritePSAModel()

registerDoParallel(cores = 4)
fixedSplines = list(c(1,2,4), c(0.5, 1.2, 2.5), c(0.5, 1.2, 2.5, 4), 
                    c(0.1, 0.5, 1.5, 4), c(0.1, 0.5, 1.5, 5), c(0.3, 0.85, 1.6, 4),
                    c(0.1, 0.5, 4), c(0.5, 4))

fixedSplines = lapply(1:50, function(i){
  first = runif(n=1, 0, 0.4)
  second = runif(n=1, first, 1)
  c(first, second, 4)
})

psaModel_splines = foreach(i=1:8, .packages = c("splines", "nlme", "ggplot2"), 
                           .export = c("psa_data_set")) %dopar%{
  psaModel_spline_trimmed = fitUnivaritePSAModel(fixedSplineKnots = fixedSplines[[i]],
                                                 randomSplineKnots = fixedSplines[[i]][1],
                                                 boundaryKnots = c(0, 7))
}

which.min(sapply(psaModel_splines, BIC))
which.min(sapply(psaModel_splines, AIC))

plotPSAFittedCurve(psaModel_splines[order(unlist(lapply(psaModel_splines, BIC)))[c(1:2,7)]], transformPSA = T, individually=F)

#######################################################

mvglm_psa_spline_124_1 = mvglmer(list(
  logpsa1 ~  I(Age - 70) +  I((Age - 70)^2) + 
    ns(visitTimeYears, knots=c(1, 2, 4), Boundary.knots=c(0,7)) + 
    (ns(visitTimeYears, knots=c(1), Boundary.knots=c(0,7))|P_ID)), 
  data =psa_data_set, families = list(gaussian))

save.image(file = "Rdata/Gleason as event/psa_gl_spline_124.Rdata")

mvJoint_psa_spline_124_1_tdval = mvJointModelBayes(mvglm_psa_spline_124_1, coxModel, timeVar = "visitTimeYears")
save.image(file = "Rdata/Gleason as event/psa_gl_spline_124.Rdata")

mvJoint_psa_spline_124_1_tdboth = mvJointModelBayes(mvglm_psa_spline_124_1, coxModel, timeVar = "visitTimeYears", 
                                   Formulas = list("logpsa1" = "value",
                                                   "logpsa1" = list(fixed = ~ 0 + dns(visitTimeYears, knots=c(1, 2, 4), Boundary.knots=c(0, 7)),
                                                                    random=~0 + dns(visitTimeYears, knots=c(1), Boundary.knots=c(0, 7)),
                                                                    indFixed = 4:7, indRandom=2:3, name = "slope")))
save.image(file = "Rdata/Gleason as event/psa_gl_spline_124.Rdata")

###############################################################
mvglm_psa_spline_124_1b15 = mvglmer(list(
  logpsa1 ~  I(Age - 70) +  I((Age - 70)^2) + 
    ns(visitTimeYears, knots=c(1, 2, 4), Boundary.knots=c(0,15)) + 
    (ns(visitTimeYears, knots=c(1), Boundary.knots=c(0,15))|P_ID)), 
  data =psa_data_set, families = list(gaussian))

save.image(file = "Rdata/Gleason as event/psa_gl_spline_124.Rdata")

mvJoint_psa_spline_124_1b15_tdval = mvJointModelBayes(mvglm_psa_spline_124_1b15, coxModel, timeVar = "visitTimeYears")
save.image(file = "Rdata/Gleason as event/psa_gl_spline_124.Rdata")

mvJoint_psa_spline_124_1b15_tdboth = mvJointModelBayes(mvglm_psa_spline_124_1b15, coxModel, timeVar = "visitTimeYears", 
                                                    Formulas = list("logpsa1" = "value",
                                                                    "logpsa1" = list(fixed = ~ 0 + dns(visitTimeYears, knots=c(1, 2, 4), Boundary.knots=c(0, 15)),
                                                                                     random=~0 + dns(visitTimeYears, knots=c(1), Boundary.knots=c(0, 15)),
                                                                                     indFixed = 4:7, indRandom=2:3, name = "slope")))
save.image(file = "Rdata/Gleason as event/psa_gl_spline_124.Rdata")


##############################################################
psa_data_set$logepsa1 = log(psa_data_set$psa + 1)

mvglm_psa_spline_124_1b15loge = mvglmer(list(
  logepsa1 ~  I(Age - 70) +  I((Age - 70)^2) + 
    ns(visitTimeYears, knots=c(1, 2, 4), Boundary.knots=c(0,15)) + 
    (ns(visitTimeYears, knots=c(1), Boundary.knots=c(0,15))|P_ID)), 
  data =psa_data_set, families = list(gaussian))

save.image(file = "Rdata/Gleason as event/psa_gl_spline_124.Rdata")

mvJoint_psa_spline_124_1b15loge_tdval = mvJointModelBayes(mvglm_psa_spline_124_1b15loge, coxModel, timeVar = "visitTimeYears")
save.image(file = "Rdata/Gleason as event/psa_gl_spline_124.Rdata")

mvJoint_psa_spline_124_1b15loge_tdboth = mvJointModelBayes(mvglm_psa_spline_124_1b15loge, coxModel, timeVar = "visitTimeYears", 
                                                       Formulas = list("logepsa1" = "value",
                                                                       "logepsa1" = list(fixed = ~ 0 + dns(visitTimeYears, knots=c(1, 2, 4), Boundary.knots=c(0, 15)),
                                                                                        random=~0 + dns(visitTimeYears, knots=c(1), Boundary.knots=c(0, 15)),
                                                                                        indFixed = 4:7, indRandom=2:3, name = "slope")))
save.image(file = "Rdata/Gleason as event/psa_gl_spline_124.Rdata")


##############################################################

mvglm_psa_spline_pt1pt54_pt1 = mvglmer(list(
  logpsa1 ~  I(Age - 70) +  I((Age - 70)^2) + 
    ns(visitTimeYears, knots=c(0.1, 0.5, 4), Boundary.knots=c(0, 7)) + 
    (ns(visitTimeYears, knots=c(0.1), Boundary.knots=c(0, 7))|P_ID)), 
  data =psa_data_set, families = list(gaussian))

save.image(file = "Rdata/Gleason as event/psa_spline_pt1pt54_pt1.Rdata")

mvJoint_psa_spline_pt1pt54_pt1_tdval = mvJointModelBayes(mvglm_psa_spline_pt1pt54_pt1, coxModel, timeVar = "visitTimeYears")
save.image(file = "Rdata/Gleason as event/psa_spline_pt1pt54_pt1.Rdata")

mvJoint_psa_spline_pt1pt54_pt11_tdboth = mvJointModelBayes(mvglm_psa_spline_pt1pt54_pt1, coxModel, timeVar = "visitTimeYears", 
                                                    Formulas = list("logpsa1" = "value",
                                                                    "logpsa1" = list(fixed = ~ 0 + dns(visitTimeYears, knots=c(0.1, 0.5, 4), Boundary.knots=c(0, 7)),
                                                                                     random=~0 + dns(visitTimeYears, knots=c(0.1), Boundary.knots=c(0, 7)),
                                                                                     indFixed = 4:7, indRandom=2:3, name = "slope")))
save.image(file = "Rdata/Gleason as event/psa_spline_pt1pt54_pt1.Rdata")
###############################################################
mvglm_psa_spline_pt51pt22pt5_pt5 = mvglmer(list(
  logpsa1 ~  I(Age - 70) +  I((Age - 70)^2) + 
    ns(visitTimeYears, knots=c(0.5, 1.2, 2.5), Boundary.knots=c(0, 7)) + 
    (ns(visitTimeYears, knots=c(0.5), Boundary.knots=c(0, 7))|P_ID)), 
  data =psa_data_set, families = list(gaussian))

save.image(file = "Rdata/Gleason as event/psa_spline_pt51pt22pt5_pt5.Rdata")

mvJoint_psa_spline_pt51pt22pt5_pt5_tdval = mvJointModelBayes(mvglm_psa_spline_pt51pt22pt5_pt5, coxModel, timeVar = "visitTimeYears")
save.image(file = "Rdata/Gleason as event/psa_spline_pt51pt22pt5_pt5.Rdata")

mvJoint_psa_spline_pt51pt22pt5_pt5_tdboth = mvJointModelBayes(mvglm_psa_spline_pt51pt22pt5_pt5, coxModel, timeVar = "visitTimeYears", 
                                                           Formulas = list("logpsa1" = "value",
                                                                           "logpsa1" = list(fixed = ~ 0 + dns(visitTimeYears, knots=c(0.5, 1.2, 2.5), Boundary.knots=c(0, 7)),
                                                                                            random=~0 + dns(visitTimeYears, knots=c(0.5), Boundary.knots=c(0, 7)),
                                                                                            indFixed = 4:7, indRandom=2:3, name = "slope")))

joint_psa_spline_pt51pt22pt5_pt5_tdboth = jointModelBayes(psaModel_splines[[2]], coxModel, timeVar = "visitTimeYears",
                                   param = "td-both",
                                   extraForm = list(fixed = ~ 0 + dns(visitTimeYears, knots=c(0.5, 1.2, 2.5), Boundary.knots=c(0, 7)),
                                                    random=~0 + dns(visitTimeYears, knots=c(0.5), Boundary.knots=c(0, 7)),
                                                    indFixed = 4:7, indRandom=2:3), control=list(n.iter = 1000))

save.image(file = "Rdata/Gleason as event/psa_spline_pt51pt22pt5_pt5.Rdata")
###############################################################
# the data frame that contains the combination of values to
# create the plot
#Hazard ratio plot
sim_psa_slope = seq(1, 1 * 5, by=1)
simhr = exp(4.583786593 * (sim_psa_slope-sim_psa_slope[1]))

qplot(x = sim_psa_slope, y=simhr, geom = "line", ylab = "Hazard Ratio", xlab = "PSA Multiplier") + 
  ticksX(from=sim_psa_slope[1], max = max(sim_psa_slope), by = sim_psa_slope[1]*2) + ticksY(from=0, 4.5, by = 0.25)
