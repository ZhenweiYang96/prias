load("Rdata/decision_analytic/cleandata.Rdata")
source("src/decision_analytic/load_lib.R")

dre_psa_data_set =  prias_long[!(is.na(prias_long$dre) & is.na(prias_long$psa)),]
dre_psa_data_set$high_dre = ifelse(dre_psa_data_set$dre=="T1c", 0, 1)

length(unique(droplevels(dre_psa_data_set$P_ID))) == length(prias.id$P_ID)

startTime_mvglmer = Sys.time()
mvglmer_dre_psa = mvglmer(list(high_dre~I(Age - 70) +  I((Age - 70)^2) + visitTimeYears  + 
                                 (visitTimeYears|P_ID),
                               
                               log2psaplus1 ~ I(Age - 70) +  I((Age - 70)^2) +
                                 ns(visitTimeYears, knots=c(0.1, 0.7, 4), Boundary.knots=c(0, 5.42)) + 
                                 (ns(visitTimeYears, knots=c(0.1, 0.7, 4), Boundary.knots=c(0, 5.42))|P_ID)),
                          
                          data=dre_psa_data_set, families = list(binomial, gaussian), engine = "STAN",
                          control = list(n.iter=10000))
endTime_mvglmer = Sys.time()

save(mvglmer_dre_psa, file="Rdata/decision_analytic/DRE_PSA/mvglmer_dre_psa.Rdata")

survModel = survreg(Surv(progression_time_start, progression_time_end, type = "interval2") ~ 
                      I(Age - 70) +  I((Age - 70)^2), data = prias.id, model = TRUE)
save(survModel, file="Rdata/decision_analytic/DRE_PSA/survModel.Rdata")

forms_dre_value = list("high_dre" = "value",
                       "log2psaplus1" = "value",
                       "log2psaplus1" = list(fixed = ~ 0 + dns(visitTimeYears, knots=c(0.1, 0.7, 4), Boundary.knots=c(0, 5.42)),
                                             random=~0 + dns(visitTimeYears, knots=c(0.1, 0.7, 4), Boundary.knots=c(0, 5.42)),
                                             indFixed = 4:7, indRandom=2:5, name = "slope"))

forms_dre_slope = list("high_dre" = "value",
                       "high_dre" = list(fixed = ~ 0 + 1, random=~0 + 1, indFixed = 4, indRandom=2, name = "slope"),
                       "log2psaplus1" = "value",
                       "log2psaplus1" = list(fixed = ~ 0 + dns(visitTimeYears, knots=c(0.1, 0.7, 4), Boundary.knots=c(0, 5.42)),
                                             random=~0 + dns(visitTimeYears, knots=c(0.1, 0.7, 4), Boundary.knots=c(0, 5.42)),
                                             indFixed = 4:7, indRandom=2:5, name = "slope"))

startTime_mvJoint_dre_value = Sys.time()
mvJoint_dre_psa_dre_value = mvJointModelBayes(mvglmer_dre_psa, survModel, 
                                              timeVar = "visitTimeYears", Formulas = forms_dre_value)
endTime_mvJoint_dre_value = Sys.time()

save(mvJoint_dre_psa_dre_value, file="Rdata/decision_analytic/DRE_PSA/mvJoint_dre_psa_dre_value.Rdata")
rm(mvJoint_dre_psa_dre_value)

##
startTime_mvJoint_dre_slope = Sys.time()
mvJoint_dre_psa_dre_slope = mvJointModelBayes(mvglmer_dre_psa, survModel, 
                                              timeVar = "visitTimeYears", Formulas = forms_dre_slope)
endTime_mvJoint_dre_slope = Sys.time()

save(mvJoint_dre_psa_dre_slope, file="Rdata/decision_analytic/DRE_PSA/mvJoint_dre_psa_dre_slope.Rdata")

