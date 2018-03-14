source("src/R/common.R")

dre_data_set =  prias_long[!is.na(prias_long$dre),]
dre_data_set$high_dre = ifelse(dre_data_set$dre=="T1c", 0, 1)

dre_data_set$Age = dre_data_set$Age - mean(dre_data_set$Age)

ggplot(data=dre_data_set, aes(y=high_dre, x=visitTimeYears)) + 
  geom_point() + stat_smooth(geom="line", colour="blue") + ticksX(from = 0, max = 20, 0.5)+ 
  xlab("Time (years)") + ylab("Probability (DRE > T1c)")

glmer_dre = glmer(high_dre ~ Age + bs(visitTimeYears, df=2, degree = 1, Boundary.knots = c(0, 5.21)) + (1|P_ID), data = dre_data_set, family = binomial(link = "logit"), 
      control = glmerControl(),
      start = NULL, verbose = 0L, nAGQ = 10L)

glmer_dre_linear = glmer(high_dre ~ Age + visitTimeYears + (1|P_ID), data = dre_data_set, family = binomial(link = "logit"), 
                  control = glmerControl(),
                  start = NULL, verbose = 0L, nAGQ = 10L)

glmer_dre_linear2 = glmer(high_dre ~ Age + visitTimeYears + (visitTimeYears|P_ID), data = dre_data_set, family = binomial(link = "logit"), 
                         control = glmerControl(),
                         start = NULL, verbose = 0L, nAGQ = 1L)


mvglmer_dre_linear = mvglmer(list(high_dre~Age +  visitTimeYears + (visitTimeYears|P_ID)),
        data = dre_data_set, families = list(binomial), engine="JAGS")

mvglmer_dre = mvglmer(list(high_dre~Age + bs(visitTimeYears, df=2, degree = 1, Boundary.knots = c(0, 5.21)) + (1|P_ID)),
                      data = dre_data_set, families = list(binomial), engine="JAGS")


mvJoint_dre_spline_autosel_tdval = mvJointModelBayes(mvglmer_dre, survModel, timeVar = "visitTimeYears")
mvJoint_dre_spline_autosel_tdboth = mvJointModelBayes(mvglmer_dre, survModel, timeVar = "visitTimeYears",
                                                      Formulas = list("high_dre" = "value",
                "high_dre" = list(fixed = ~ 0 + dbs_deg1(visitTimeYears, df=2, Boundary.knots=c(0, 7)),
                                 random=~0,
                                 indFixed = 3:4, name = "slope")))

mvJoint_dre_linear_tdboth = mvJointModelBayes(mvglmer_dre_linear, survModel, timeVar = "visitTimeYears",
                                                      Formulas = list("high_dre" = "value",
                                                                      "high_dre" = list(fixed = ~ 0 + 1,
                                                                                        random=~0 + 1,
                                                                                        indFixed = 3, indRandom=2, name = "slope")))

