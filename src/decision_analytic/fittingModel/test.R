mvglmer_psa_simple_randeff_normaldist = mvglmer(list(log2psa ~ I(Age - 70) +  I((Age - 70)^2) + 
                                            ns(visitTimeYears, knots=c(0.1, 0.7, 4), Boundary.knots=c(0, 5.42)) + 
                                            (ns(visitTimeYears, knots=c(0.1), Boundary.knots=c(0, 5.42))|P_ID)),
                                     data=psa_data_set, families = list(gaussian))
endTime_mvglmer = Sys.time()

save(mvglmer_psa_simple_randeff_normaldist, file="Rdata/decision_analytic/PSA_Only/mvglmer_psa_simple_randeff_normaldist.Rdata")
forms = list("log2psa" = "value",
             "log2psa" = list(fixed = ~ 0 + dns(visitTimeYears, knots=c(0.1, 0.7, 4), Boundary.knots=c(0, 5.42)),
                              random=~0 + dns(visitTimeYears, knots=c(0.1), Boundary.knots=c(0, 5.42)),
                              indFixed = 4:7, indRandom=2:3, name = "slope"))

startTime_mvJoint_psa = Sys.time()
mvJoint_psa_simple_randeff_normaldist = mvJointModelBayes(mvglmer_psa_simple_randeff_normaldist, survModel, timeVar = "visitTimeYears", Formulas = forms)
endTime_mvJoint_psa = Sys.time()
save(mvJoint_psa_simple_randeff_normaldist, file="Rdata/decision_analytic/PSA_Only/mvJoint_psa_simple_randeff_normaldist.Rdata")
rm(mvglmer_psa_simple_randeff_normaldist)
rm(mvJoint_psa_simple_randeff_normaldist)

startTime_mvglmer = Sys.time()
mvglmer_psa_normaldist = mvglmer(list(log2psa ~ I(Age - 70) +  I((Age - 70)^2) + 
                             ns(visitTimeYears, knots=c(0.1, 0.7, 4), Boundary.knots=c(0, 5.42)) + 
                             (ns(visitTimeYears, knots=c(0.1, 0.7, 4), Boundary.knots=c(0, 5.42))|P_ID)),
                      data=psa_data_set, families = list(gaussian))
endTime_mvglmer = Sys.time()

save(mvglmer_psa_normaldist, file="Rdata/decision_analytic/PSA_Only/mvglmer_psa_normaldist.Rdata")

survModel = survreg(Surv(progression_time_start, progression_time_end, type = "interval2") ~ 
                      I(Age - 70) +  I((Age - 70)^2), data = prias.id, model = TRUE)
save(survModel, file="Rdata/decision_analytic/PSA_Only/survModel.Rdata")

forms = list("log2psa" = "value",
             "log2psa" = list(fixed = ~ 0 + dns(visitTimeYears, knots=c(0.1, 0.7, 4), Boundary.knots=c(0, 5.42)),
                              random=~0 + dns(visitTimeYears, knots=c(0.1, 0.7, 4), Boundary.knots=c(0, 5.42)),
                              indFixed = 4:7, indRandom=2:5, name = "slope"))

startTime_mvJoint_psa = Sys.time()
mvJoint_psa_normaldist = mvJointModelBayes(mvglmer_psa_normaldist, survModel, timeVar = "visitTimeYears", Formulas = forms)
endTime_mvJoint_psa = Sys.time()
save(mvJoint_psa_normaldist, file="Rdata/decision_analytic/PSA_Only/mvJoint_psa_normaldist.Rdata")
rm(mvglmer_psa_normaldist)
rm(mvJoint_psa_normaldist)

load("Rdata/decision_analytic/DRE Only/mvglmer_dre_linear_random_slope_full.Rdata")
mvJoint_dre_linear_tdval = mvJointModelBayes(mvglmer_dre_linear_random_slope_full, survModel, timeVar = "visitTimeYears")
save(mvJoint_dre_linear_tdval, file="Rdata/decision_analytic/DRE Only/mvJoint_dre_linear_tdval.Rdata")


mvJoint_dre_linear_tdboth = mvJointModelBayes(mvglmer_dre_linear_random_slope_full, survModel, timeVar = "visitTimeYears",
                                              Formulas = list("high_dre" = "value",
                                                              "high_dre" = list(fixed = ~ 0 + 1,
                                                                                random=~0 + 1,
                                                                                indFixed = 4, indRandom=2, name = "slope")))
save(mvJoint_dre_linear_tdboth, file="Rdata/decision_analytic/DRE Only/mvJoint_dre_linear_tdboth.Rdata")



###########################
survModel = survreg(Surv(progression_time_start, progression_time_end, type = "interval2") ~ 
                      I(Age - 70) +  I((Age - 70)^2), data = prias.id, model = TRUE)
save(survModel, file="Rdata/decision_analytic/DRE_PSA/survModel.Rdata")

forms_dre_value = list("high_dre" = "value",
                                    "log2psa" = "value",
                                    "log2psa" = list(fixed = ~ 0 + dns(visitTimeYears, knots=c(0.1, 0.7, 4), Boundary.knots=c(0, 5.42)),
                                                     random=~0 + dns(visitTimeYears, knots=c(0.1, 0.7, 4), Boundary.knots=c(0, 5.42)),
                                                     indFixed = 4:7, indRandom=2:5, name = "slope"))

forms_dre_slope = list("high_dre" = "value",
                                    "high_dre" = list(fixed = ~ 0 + 1, random=~0 + 1, indFixed = 4, indRandom=2, name = "slope"),
                                    "log2psa" = "value",
                                    "log2psa" = list(fixed = ~ 0 + dns(visitTimeYears, knots=c(0.1, 0.7, 4), Boundary.knots=c(0, 5.42)),
                                                     random=~0 + dns(visitTimeYears, knots=c(0.1, 0.7, 4), Boundary.knots=c(0, 5.42)),
                                                     indFixed = 4:7, indRandom=2:5, name = "slope"))

startTime_mvJoint_dre_value = Sys.time()
mvJoint_dre_psa_dre_value = mvJointModelBayes(mvglmer_dre_psa, survModel, 
                                                           timeVar = "visitTimeYears", Formulas = forms_dre_value)
endTime_mvJoint_dre_value = Sys.time()

save(mvJoint_dre_psa_dre_value, file="Rdata/decision_analytic/DRE_PSA/mvJoint_dre_psa_dre_value.Rdata")
rm(mvJoint_dre_psa_dre_value)

###########
startTime_mvJoint_dre_value_expit = Sys.time()
mvJoint_dre_psa_dre_value_expit = mvJointModelBayes(mvglmer_dre_psa, survModel, 
                                                                 timeVar = "visitTimeYears", Formulas = forms_dre_value,
                                                                 transFuns = c("high_dre_value"="expit", "log2psa_value"="identity", "log2psa_slope"="identity"))
endTime_mvJoint_dre_value_expit = Sys.time()

save(mvJoint_dre_psa_dre_value_expit, file="Rdata/decision_analytic/DRE_PSA/mvJoint_dre_psa_dre_value_expit.Rdata")
rm(mvJoint_dre_psa_dre_value_expit)

startTime_mvJoint_dre_slope = Sys.time()
mvJoint_dre_psa_dre_slope = mvJointModelBayes(mvglmer_dre_psa, survModel, 
                                              timeVar = "visitTimeYears", Formulas = forms_dre_slope)
endTime_mvJoint_dre_slope = Sys.time()

save(mvJoint_dre_psa_dre_slope, file="Rdata/decision_analytic/DRE_PSA/mvJoint_dre_psa_dre_slope.Rdata")
rm(mvJoint_dre_psa_dre_slope)