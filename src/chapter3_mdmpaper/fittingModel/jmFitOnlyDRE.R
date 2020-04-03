load("Rdata/decision_analytic/cleandata.Rdata")
source("src/decision_analytic/load_lib.R")

dre_data_set =  prias_long[!is.na(prias_long$dre),]
dre_data_set$high_dre = ifelse(dre_data_set$dre=="T1c", 0, 1)

ggplot(data=dre_data_set, aes(y=high_dre, x=visitTimeYears)) + 
  geom_point() + stat_smooth(geom="line", colour="blue") + 
  xlab("Time (years)") + ylab("Probability (DRE > T1c)")

dre_data_set$std_age = (dre_data_set$Age - 70)/sd(dre_data_set$Age)
dre_data_set$std_visitTimeYears = (dre_data_set$visitTimeYears - mean(dre_data_set$visitTimeYears))/(sd(dre_data_set$visitTimeYears))
glmer_dre = glmer(high_dre ~ std_age +  I(std_age^2) + 
                    std_visitTimeYears  + (std_visitTimeYears|P_ID), 
                  data = dre_data_set, family = binomial(link = "logit"), 
                         control = glmerControl())

mvglmer_dre_linear_random_intercept = mvglmer(list(high_dre~I(Age - 70) +  I((Age - 70)^2) + 
                                    visitTimeYears + (1|P_ID)),
        data = dre_data_set, families = list(binomial), engine="JAGS")
save(mvglmer_dre_linear_random_intercept, file = "Rdata/decision_analytic/mvglmer_dre_linear_random_intercept.Rdata")

mvglmer_dre = mvglmer(list(high_dre~ I(Age - 70) + I((Age - 70)^2) + visitTimeYears + (visitTimeYears|P_ID)),
                                               data = dre_data_set, families = list(binomial), engine="JAGS")
save(mvglmer_dre, file = "Rdata/decision_analytic/DRE Only/mvglmer_dre.Rdata")

survModel = survreg(Surv(progression_time_start, progression_time_end, type = "interval2") ~ 
                      I(Age - 70) +  I((Age - 70)^2), data = prias.id, model = TRUE)
save(survModel, file="Rdata/decision_analytic/DRE Only/survModel.Rdata")

mvJoint_dre_tdval = mvJointModelBayes(mvglmer_dre, survModel, timeVar = "visitTimeYears")
save(mvJoint_dre_tdval, file="Rdata/decision_analytic/DRE Only/mvJoint_dre_tdval.Rdata")

mvJoint_dre_tdboth = mvJointModelBayes(mvglmer_dre, survModel, timeVar = "visitTimeYears",
                                                      Formulas = list("high_dre" = "value",
                                                                      "high_dre" = list(fixed = ~ 0 + 1,
                                                                                        random=~0 + 1,
                                                                                       indFixed = 4, indRandom=2, name = "slope")))
save(mvJoint_dre_tdboth, file="Rdata/decision_analytic/DRE Only/mvJoint_dre_tdboth.Rdata")

