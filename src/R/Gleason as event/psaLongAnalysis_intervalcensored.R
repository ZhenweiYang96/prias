source("src/R/common.R")

psa_data_set =  prias_long[!is.na(prias_long$log2psa),]

pidList = unique(psa_data_set$P_ID)
randomSample = sample(x = pidList, replace = F, size = 150)
ggplot(data=psa_data_set, aes(x=visitTimeYears, y=log2psa)) + 
  geom_line(aes(group=P_ID)) + stat_smooth() + 
  ticksX(from=0,max(psa_data_set$visitTimeYears), 0.25) + xlab("Time (years)") + 
  ylab("log(psa): Adjusted for Age")


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

########################################################################
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
