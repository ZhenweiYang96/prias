source("src/R/common.R")

training_long = prias_long[!prias_long$P_ID %in% c(3174, 2340, 911),]

psa_data_set =  prias_long[!is.na(prias_long$log2psa),]
training_psa_data_set = training_long[!is.na(training_long$log2psa),]

training_psa_data_set$progression_time = ifelse(training_psa_data_set$progression_time_end==Inf, 
                                                training_psa_data_set$progression_time_start, 
                                                training_psa_data_set$progression_time_end)


ggplot(data=psa_data_set[psa_data_set$P_ID==3174,], aes(x=visitTimeYears, y=psa)) + 
  geom_line() + geom_point() + 
  ticksX(from=0,max(psa_data_set$visitTimeYears), 1) + xlab("Time (years)") + 
  ylab("PSA (ng/mL)") + ggtitle(paste("Patient ID:", 3174))

########################################################################
ct= makeCluster(detectCores())
registerDoParallel(ct)

fixedSplines = list(c(1,2,4), c(0.5, 1.2, 2.5), c(0.5, 1.2, 2.5, 4), 
                    c(0.1, 0.5, 1.5, 4), c(0.1, 0.5, 1.5, 5), c(0.3, 0.85, 1.6, 4),
                    c(0.1, 0.5, 4), c(0.5, 4))

fixedSplines = lapply(1:50, function(i){
  first = runif(n=1, 0, 0.4)
  second = runif(n=1, first, 1)
  c(first, second, 4)
})

psaModel_splines = foreach(i=1:length(fixedSplines), .packages = c("splines", "nlme", "ggplot2"), 
                           .export = c("psa_data_set")) %dopar%{
                             psaModel_spline_trimmed = fitUnivaritePSAModel(fixedSplineKnots = fixedSplines[[i]],
                                                                            randomSplineKnots = fixedSplines[[i]][1],
                                                                            boundaryKnots = c(0, 7))
                           }


stopCluster(ct)


which.min(sapply(psaModel_splines, BIC))
which.min(sapply(psaModel_splines, AIC))

plotPSAFittedCurve(psaModel_splines[c(1:2,7)], transformPSA = T, individually=F)

ggplot(data=psa_data_set[psa_data_set$visitTimeYears>4,], aes(x=visitTimeYears, y=log2psa)) + 
  geom_line(aes(group=P_ID)) + stat_smooth() + 
  ticksX(from=0,max(psa_data_set$visitTimeYears), 1) + xlab("Time (years)") + 
  ylab(expression('log'[2]*'(PSA)'))

########################################################################
mvglm_psa_spline_pt1pt54_pt1 = mvglmer(list(
  log2psa ~  I(Age - 70) +  I((Age - 70)^2) + 
    ns(visitTimeYears, knots=c(0.1, 0.5, 4), Boundary.knots=c(0, 7)) + 
    (ns(visitTimeYears, knots=c(0.1), Boundary.knots=c(0, 7))|P_ID)), 
  data=training_psa_data_set, families = list(gaussian))

save.image(file = "Rdata/Gleason as event/psa_spline_pt1pt54_pt1.Rdata")

mvJoint_psa_spline_pt51pt22pt5_pt5_tdval = mvJointModelBayes(mvglm_psa_spline_pt1pt54_pt1, survModel.training, 
                                                             timeVar = "visitTimeYears")
save.image(file = "Rdata/Gleason as event/psa_spline_pt1pt54_pt1.Rdata")

mvJoint_psa_spline_pt1pt54_pt1_tdboth = mvJointModelBayes(mvglm_psa_spline_pt1pt54_pt1, survModel.training, timeVar = "visitTimeYears", 
                                                              Formulas = list("log2psa" = "value",
                                                                              "log2psa" = list(fixed = ~ 0 + dns(visitTimeYears, knots=c(0.1, 0.5, 4), Boundary.knots=c(0, 7)),
                                                                                               random=~0 + dns(visitTimeYears, knots=c(0.1), Boundary.knots=c(0, 7)),
                                                                                               indFixed = 4:7, indRandom=2:3, name = "slope")))

lme_psa_spline_pt1pt54_pt1 = lme(fixed=log2psa ~  I(Age - 70) +  I((Age - 70)^2) + 
                                   ns(visitTimeYears, knots=c(0.1, 0.5, 4), Boundary.knots=c(0, 7)), 
                                 random = ~ns(visitTimeYears, knots=c(0.1), Boundary.knots=c(0, 7))|P_ID, 
                                 data=training_psa_data_set,
                                 control = lmeControl(opt = "optim", optimMethod = "L-BFGS-B"), 
                                 method = "REML")

joint_psa_spline_pt1pt54_pt1_tdboth = jointModelBayes(lme_psa_spline_pt1pt54_pt1, survModel.training_rightCens, 
                                                      timeVar = "visitTimeYears", param = "td-both",
                                                      extraForm = list(fixed = ~ 0 + dns(visitTimeYears, knots=c(0.1, 0.5, 4), Boundary.knots=c(0, 7)),
                                                                       random=~0 + dns(visitTimeYears, knots=c(0.1), Boundary.knots=c(0, 7)),
                                                                       indFixed = 4:7, indRandom=2:3), control=list(n.iter = 1000))

save.image(file = "Rdata/Gleason as event/psa_spline_pt1pt54_pt1.Rdata")

joint_psa_replaced_prias = replaceMCMCContents(mvJoint_psa_spline_pt1pt54_pt1_tdboth, 
                                               joint_psa_spline_pt1pt54_pt1_tdboth)
save.image(file = "Rdata/Gleason as event/psa_spline_pt1pt54_pt1.Rdata")

###############for right censoring data#################
mvJoint_psa_spline_pt1pt54_pt1_tdboth_rc = mvJointModelBayes(mvglm_psa_spline_pt1pt54_pt1, survModel.training_rightCens, timeVar = "visitTimeYears", 
                                                             Formulas = list("log2psa" = "value",
                                                                             "log2psa" = list(fixed = ~ 0 + dns(visitTimeYears, knots=c(0.1, 0.5, 4), Boundary.knots=c(0, 7)),
                                                                                              random=~0 + dns(visitTimeYears, knots=c(0.1), Boundary.knots=c(0, 7)),
                                                                                              indFixed = 4:7, indRandom=2:3, name = "slope")))

save.image(file = "Rdata/Gleason as event/psa_spline_pt1pt54_pt1_rc.Rdata")

joint_psa_replaced_prias_rc = replaceMCMCContents(mvJoint_psa_spline_pt1pt54_pt1_tdboth_rc, 
                                               joint_psa_spline_pt1pt54_pt1_tdboth)
save.image(file = "Rdata/Gleason as event/psa_spline_pt1pt54_pt1_rc.Rdata")

########################################################################
mvglm_psa_spline_pt51pt22pt5_pt5 = mvglmer(list(
  log2psa ~  I(Age - 70) +  I((Age - 70)^2) + 
    ns(visitTimeYears, knots=c(0.5, 1.2, 2.5), Boundary.knots=c(0, 7)) + 
    (ns(visitTimeYears, knots=c(0.5), Boundary.knots=c(0, 7))|P_ID)), 
  data =psa_data_set, families = list(gaussian))

save.image(file = "Rdata/Gleason as event/psa_spline_pt51pt22pt5_pt5.Rdata")

mvJoint_psa_spline_pt51pt22pt5_pt5_tdval = mvJointModelBayes(mvglm_psa_spline_pt51pt22pt5_pt5, survModel, timeVar = "visitTimeYears")
save.image(file = "Rdata/Gleason as event/psa_spline_pt51pt22pt5_pt5.Rdata")

mvJoint_psa_spline_pt51pt22pt5_pt5_tdboth = mvJointModelBayes(mvglm_psa_spline_pt51pt22pt5_pt5, survModel, timeVar = "visitTimeYears", 
                                                           Formulas = list("log2psa" = "value",
                                                                           "log2psa" = list(fixed = ~ 0 + dns(visitTimeYears, knots=c(0.5, 1.2, 2.5), Boundary.knots=c(0, 7)),
                                                                                            random=~0 + dns(visitTimeYears, knots=c(0.5), Boundary.knots=c(0, 7)),
                                                                                            indFixed = 4:7, indRandom=2:3, name = "slope")))

lme_psa_spline_pt51pt22pt5_pt5 = lme(fixed=log2psa ~  I(Age - 70) +  I((Age - 70)^2) + 
              ns(visitTimeYears, knots=c(0.5, 1.2, 2.5), Boundary.knots=c(0, 7)), 
            random = ~ns(visitTimeYears, knots=c(0.5), Boundary.knots=c(0, 7))|P_ID, 
            data=psa_data_set,
            control = lmeControl(opt = "optim", optimMethod = "L-BFGS-B"), 
            method = "REML")

joint_psa_spline_pt51pt22pt5_pt5_tdboth = jointModelBayes(lme_psa_spline_pt51pt22pt5_pt5, survModel_rightCens, 
                                                          timeVar = "visitTimeYears", param = "td-both",
                                   extraForm = list(fixed = ~ 0 + dns(visitTimeYears, knots=c(0.5, 1.2, 2.5), Boundary.knots=c(0, 7)),
                                                    random=~0 + dns(visitTimeYears, knots=c(0.5), Boundary.knots=c(0, 7)),
                                                    indFixed = 4:7, indRandom=2:3), control=list(n.iter = 1000))

save.image(file = "Rdata/Gleason as event/psa_spline_pt51pt22pt5_pt5.Rdata")
########################################################################
joint_psa_replaced_prias = replaceMCMCContents(mvJoint_psa_spline_pt51pt22pt5_pt5_tdboth, 
                                               joint_psa_spline_pt51pt22pt5_pt5_tdboth)
save.image(file = "Rdata/Gleason as event/psa_spline_pt51pt22pt5_pt5.Rdata")

#########################################################################
