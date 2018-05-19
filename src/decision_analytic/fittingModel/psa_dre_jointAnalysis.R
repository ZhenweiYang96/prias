source("src/decision_analytic/load_lib.R")

dre_psa_data_set =  prias_long[!(is.na(prias_long$dre) & is.na(prias_long$psa)),]
dre_psa_data_set$high_dre = ifelse(dre_psa_data_set$dre=="T1c", 0, 1)

length(unique(droplevels(dre_psa_data_set$P_ID))) == length(prias.id$P_ID)

#Standardize the variables
dre_psa_data_set$age_std = scale(dre_psa_data_set$Age, center = T, scale = T)
prias.id$age_std = dre_psa_data_set$age_std[!duplicated(dre_psa_data_set$P_ID)]

mean(dre_psa_data_set$visitTimeYears)
sd(dre_psa_data_set$visitTimeYears)
std_knots = (c(0.1, 0.7, 4) - mean(dre_psa_data_set$visitTimeYears))/sd(dre_psa_data_set$visitTimeYears)
std_boundary_knots = (c(0, 5.42) - mean(dre_psa_data_set$visitTimeYears))/sd(dre_psa_data_set$visitTimeYears)

##########CODE RUNNING OVER A DAY ###############
source("src/decision_analytic/load_t_dist_JMbayes.R")
startTime_mvglmer = Sys.time()
mvglmer_dre_psa_standardized = mvglmer(list(high_dre~age_std + I(age_std^2) +
                                          I((visitTimeYears - 1.812641)/1.727464) + 
                                            (I((visitTimeYears - 1.812641)/1.727464)|P_ID),
                                        
                                        log2psa ~ age_std + I(age_std^2) +
                                          ns(I((visitTimeYears - 1.812641)/1.727464), knots=c(-0.9914195, -0.6440895,  1.2662257), Boundary.knots=c(-1.049308, 2.088240)) + 
                                          (ns(I((visitTimeYears - 1.812641)/1.727464), knots=c(-0.9914195, -0.6440895,  1.2662257), Boundary.knots=c(-1.049308, 2.088240))|P_ID)),
                                   
                                   data=dre_psa_data_set, families = list(binomial, gaussian), engine = "STAN",
                                   control = list(n.iter=2000))
endTime_mvglmer = Sys.time()

save(mvglmer_dre_psa_standardized, file="Rdata/decision_analytic/DRE_PSA/mvglmer_dre_psa_standardized.Rdata")

startTime_mvglmer = Sys.time()
mvglmer_dre_psa = mvglmer(list(high_dre~I(Age - 70) +  I((Age - 70)^2) + visitTimeYears  + 
                                              (visitTimeYears|P_ID),
                                            
                                            log2psa ~ I(Age - 70) +  I((Age - 70)^2) +
                                              ns(visitTimeYears, knots=c(0.1, 0.7, 4), Boundary.knots=c(0, 5.42)) + 
                                              (ns(visitTimeYears, knots=c(0.1, 0.7, 4), Boundary.knots=c(0, 5.42))|P_ID)),
                                       
                                       data=dre_psa_data_set, families = list(binomial, gaussian), engine = "STAN",
                                       control = list(n.iter=2500))
endTime_mvglmer = Sys.time()

save(mvglmer_dre_psa, file="Rdata/decision_analytic/DRE_PSA/mvglmer_dre_psa.Rdata")



survModel_standardized = survreg(Surv(progression_time_start, progression_time_end, type = "interval2") ~ 
                      age_std + I(age_std^2), data = prias.id, model = TRUE)
save(survModel_standardized, file="Rdata/decision_analytic/DRE_PSA/survModel_standardized.Rdata")

forms_dre_value_standardized = list("high_dre" = "value",
                       "log2psa" = "value",
                       "log2psa" = list(fixed = ~ 0 + dns(I((visitTimeYears - 1.812641)/1.727464), knots=c(-0.9914195, -0.6440895,  1.2662257), Boundary.knots=c(-1.049308, 2.088240)),
                                        random=~0 + dns(I((visitTimeYears - 1.812641)/1.727464), knots=c(-0.9914195, -0.6440895,  1.2662257), Boundary.knots=c(-1.049308, 2.088240)),
                                        indFixed = 4:7, indRandom=2:5, name = "slope"))

forms_dre_slope_standardized = list("high_dre" = "value",
                 "high_dre" = list(fixed = ~ 0 + 1, random=~0 + 1, indFixed = 4, indRandom=2, name = "slope"),
                 "log2psa" = "value",
                 "log2psa" = list(fixed = ~ 0 + dns(I((visitTimeYears - 1.812641)/1.727464), knots=c(-0.9914195, -0.6440895,  1.2662257), Boundary.knots=c(-1.049308, 2.088240)),
                                  random=~0 + dns(I((visitTimeYears - 1.812641)/1.727464), knots=c(-0.9914195, -0.6440895,  1.2662257), Boundary.knots=c(-1.049308, 2.088240)),
                                  indFixed = 4:7, indRandom=2:5, name = "slope"))

startTime_mvJoint_dre_value = Sys.time()
mvJoint_dre_psa_dre_value_standardized = mvJointModelBayes(mvglmer_dre_psa_standardized, survModel_standardized, 
                                              timeVar = "visitTimeYears", Formulas = forms_dre_value_standardized)
endTime_mvJoint_dre_value = Sys.time()

save(mvJoint_dre_psa_dre_value_standardized, file="Rdata/decision_analytic/mvJoint_dre_psa_dre_value_standardized.Rdata")
rm(mvJoint_dre_psa_dre_value_standardized)

###########
startTime_mvJoint_dre_value_expit_standardized = Sys.time()
mvJoint_dre_psa_dre_value_expit_standardized = mvJointModelBayes(mvglmer_dre_psa_standardized, survModel_standardized, 
                                              timeVar = "visitTimeYears", Formulas = forms_dre_value_standardized,
                                              transFuns = c("high_dre_value"="expit", "log2psa_value"="identity", "log2psa_slope"="identity"))
endTime_mvJoint_dre_value_expit = Sys.time()

save(mvJoint_dre_psa_dre_value_expit_standardized, file="Rdata/decision_analytic/mvJoint_dre_psa_dre_value_expit_standardized.Rdata")
rm(mvJoint_dre_psa_dre_value_expit)

startTime_mvJoint_dre_slope = Sys.time()
mvJoint_dre_psa_dre_slope_standardized = mvJointModelBayes(mvglmer_dre_psa_standardized, survModel_standardized, 
                                              timeVar = "visitTimeYears", Formulas = forms_dre_slope_standardized)
endTime_mvJoint_dre_slope = Sys.time()

save(mvJoint_dre_psa_dre_slope_standardized, file="Rdata/decision_analytic/mvJoint_dre_psa_dre_slope_standardized.Rdata")
