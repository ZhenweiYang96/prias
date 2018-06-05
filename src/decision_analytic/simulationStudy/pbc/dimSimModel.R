lmeFit <- lme(y ~ year, random = ~year|id, data = DF)
CoxFit1_training <- coxph(Surv(Time, event) ~ drug + age, data = DF.id, model = TRUE,x = T)
jmFit <- jointModel(lmeFit, CoxFit1_training, timeVar = "year", method = "weibull-PH-aGH")
jmFit <- jointModel(lmeFit, CoxFit1_training, timeVar = "year", method = "piecewise-PH-GH")

summary(jmFit)

################
jmBayesFit = jointModelBayes(lmeObject = lmeFit, survObject = CoxFit1_training, timeVar = "year")
##############
mvglmer_Fit = mvglmer(list(high_dre~I(Age - 70) +  I((Age - 70)^2) + visitTimeYears  + 
                                       (visitTimeYears|P_ID),
                                     
                                     log2psaplus1 ~ I(Age - 70) +  I((Age - 70)^2) +
                                       ns(visitTimeYears, knots=c(0.1, 0.7, 4), Boundary.knots=c(0, 5.42)) + 
                                       (ns(visitTimeYears, knots=c(0.1, 0.7, 4), Boundary.knots=c(0, 5.42))|P_ID)),
                                
                                data=trainingDs, families = list(binomial, gaussian), engine = "STAN",
                                control = list(n.iter=mvglmer_iter, n.processors=max_cores))

mvJoint_dre_psa_simDs = mvJointModelBayes(mvglmer_Fit, coxFit, timeVar = "time")
#################################################

library(JM)

lmeFit <- lme(log2psaplus1 ~ I(Age - 70) +  I((Age - 70)^2) + ns(visitTimeYears, knots=c(0.1, 0.7, 4), Boundary.knots=c(0, 5.42)),
              random = ~ns(visitTimeYears, knots=c(0.1), Boundary.knots=c(0, 5.42))|P_ID, 
              data = jointModelData$trainingData$trainingDs,
              control = lmeControl(opt = "optim", optimMethod = "L-BFGS-B"))

jmFit <- jointModel(lmeFit, jointModelData$survModel_simDs, timeVar = "visitTimeYears", 
                    method = "weibull-PH-aGH", parameterization = "both",
                    derivForm = list(fixed = ~ 0 + dns(visitTimeYears, knots=c(0.1, 0.7, 4), Boundary.knots=c(0, 5.42)),
                           random=~0 + dns(visitTimeYears, knots=c(0.1), Boundary.knots=c(0, 5.42)),
                           indFixed = 4:7, indRandom=2:3, name = "slope"))


res = list(lmeFit, jmFit)
save.image("Rdata/res.Rdata")
