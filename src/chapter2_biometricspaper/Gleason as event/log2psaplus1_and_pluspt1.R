training_psa_data_set$log2psa_plus1 = log(training_psa_data_set$psa + 1, base = 2)
training_psa_data_set$log2psa_pluspt1 = log(training_psa_data_set$psa + 0.1, base = 2)

mvglm_log2psa_plus1_spline_pt1pt54_pt1 = mvglmer(list(
  log2psa_plus1 ~ I(Age - 70) +  I((Age - 70)^2) + 
    ns(visitTimeYears, knots=c(0.1, 0.5, 4), Boundary.knots=c(0, 7)) + 
    (ns(visitTimeYears, knots=c(0.1), Boundary.knots=c(0, 7))|P_ID)), 
  data=training_psa_data_set, families = list(gaussian))

save(mvglm_log2psa_plus1_spline_pt1pt54_pt1, file = "Rdata/Gleason as event/tdist/log2psa_plus1_andplus_pt1/mvglm_log2psa_plus1_spline_pt1pt54_pt1.Rdata")

mvjoint_log2psa_plus1_spline_pt1pt54_pt1 = mvJointModelBayes(mvglm_log2psa_plus1_spline_pt1pt54_pt1, survModel.training, timeVar = "visitTimeYears", 
                                                             Formulas = list("log2psa_plus1" = "value",
                                                                             "log2psa_plus1" = list(fixed = ~ 0 + dns(visitTimeYears, knots=c(0.1, 0.5, 4), Boundary.knots=c(0, 7)),
                                                                                              random=~0 + dns(visitTimeYears, knots=c(0.1), Boundary.knots=c(0, 7)),
                                                                                              indFixed = 4:7, indRandom=2:3, name = "slope")))

save(mvjoint_log2psa_plus1_spline_pt1pt54_pt1, file = "Rdata/Gleason as event/tdist/log2psa_plus1_andplus_pt1/mvjoint_log2psa_plus1_spline_pt1pt54_pt1.Rdata")

rm(mvglm_log2psa_plus1_spline_pt1pt54_pt1)
rm(mvjoint_log2psa_plus1_spline_pt1pt54_pt1)

##########################
mvglm_log2psa_pluspt1_spline_pt1pt54_pt1 = mvglmer(list(
  log2psa_pluspt1 ~  I(Age - 70) +  I((Age - 70)^2) + 
    ns(visitTimeYears, knots=c(0.1, 0.5, 4), Boundary.knots=c(0, 7)) + 
    (ns(visitTimeYears, knots=c(0.1), Boundary.knots=c(0, 7))|P_ID)), 
  data=training_psa_data_set, families = list(gaussian))

save(mvglm_log2psa_pluspt1_spline_pt1pt54_pt1, file = "Rdata/Gleason as event/tdist/log2psa_plus1_andplus_pt1/mvglm_log2psa_pluspt1_spline_pt1pt54_pt1.Rdata")

mvjoint_log2psa_pluspt1_spline_pt1pt54_pt1 = mvJointModelBayes(mvglm_log2psa_pluspt1_spline_pt1pt54_pt1, survModel.training, timeVar = "visitTimeYears", 
                                                             Formulas = list("log2psa_pluspt1" = "value",
                                                                             "log2psa_pluspt1" = list(fixed = ~ 0 + dns(visitTimeYears, knots=c(0.1, 0.5, 4), Boundary.knots=c(0, 7)),
                                                                                              random=~0 + dns(visitTimeYears, knots=c(0.1), Boundary.knots=c(0, 7)),
                                                                                              indFixed = 4:7, indRandom=2:3, name = "slope")))

save(mvjoint_log2psa_pluspt1_spline_pt1pt54_pt1, file = "Rdata/Gleason as event/tdist/log2psa_plus1_andplus_pt1/mvjoint_log2psa_pluspt1_spline_pt1pt54_pt1.Rdata")

rm(mvglm_log2psa_pluspt1_spline_pt1pt54_pt1)
rm(mvjoint_log2psa_pluspt1_spline_pt1pt54_pt1)