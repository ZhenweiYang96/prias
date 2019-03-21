load("Rdata/decision_analytic/cleandata.Rdata")
source("src/decision_analytic/load_lib.R")

psa_data_set = prias_long[!is.na(prias_long$psa),]

ggplot(data=psa_data_set[psa_data_set$visitTimeYears <= 4,], aes(x=visitTimeYears, y=log2psa)) + 
  geom_line(aes(group=P_ID)) + stat_smooth() + 
  xlab("Time (years)") +  ylab(expression('log'[2]*'(PSA)'))

#Testing various models
firstKnot = seq(0.1, 1, by = 0.1)
secondKnot = seq(0.5, 2.5, by = 0.2)
thirdKnot = seq(1.5, 4, by=0.25)

totalModels = 0
knotsList = list()
for(i in 1:length(firstKnot)){
  for(j in 1:length(secondKnot)){
    for(k in 1:length(thirdKnot)){
      if(all(diff(c(firstKnot[i], secondKnot[j], thirdKnot[k])) > 0)){
        totalModels = totalModels + 1
        knotsList[[totalModels]] = c(firstKnot[i], secondKnot[j], thirdKnot[k])
      }
    }
  }
}

ct = makeCluster(detectCores(), type="FORK")
registerDoParallel(ct)

startTime = Sys.time()

aic_bic_list = foreach(knots = knotsList, .packages = c("lme4", "splines"), .combine = "rbind") %dopar%{
  model = lmer(log2psa ~  I(Age - 70) +  I((Age - 70)^2) + 
         ns(visitTimeYears, knots=knots, Boundary.knots=c(0, 5.42)) +
         (ns(visitTimeYears, knots=knots, Boundary.knots=c(0, 5.42))|P_ID), 
       REML=F, data = psa_data_set)
  
  return(c(AIC(model), BIC(model)))
}

endTime = Sys.time()
stopCluster(ct)

knotsList[[which.min(aic_bic_list[,1])]]
knotsList[[which.min(aic_bic_list[,2])]]
#the model with the least AIC as well as least BIC is the one with knots 0.1, 0.7, 4

psaModelTestResult = list(aic_bic_list=aic_bic_list, knotsList=knotsList)
save(psaModelTestResult, file = "Rdata/decision_analytic/PSA_Only/psaModelTestResult.Rdata")

lme(fixed=log2psa ~  I(Age - 70) +  I((Age - 70)^2) + 
                             ns(visitTimeYears, knots=c(0.1, 0.7, 4), Boundary.knots=c(0, 5.42)), 
                           random = ~ns(visitTimeYears, knots=c(0.1, 0.7, 4), Boundary.knots=c(0, 5.42))|P_ID, 
                           data=psa_data_set,
                           control = lmeControl(opt = "optim", optimMethod = "L-BFGS-B"), 
                           method = "ML")

startTime_mvglmer = Sys.time()
mvglmer_psa = mvglmer(list(log2psaplus1 ~ I(Age - 70) +  I((Age - 70)^2) + 
                             ns(visitTimeYears, knots=c(0.1, 0.7, 4), Boundary.knots=c(0, 5.42)) + 
                                 (ns(visitTimeYears, knots=c(0.1, 0.7, 4), Boundary.knots=c(0, 5.42))|P_ID)),
                          data=psa_data_set, families = list(gaussian))
endTime_mvglmer = Sys.time()

save(mvglmer_psa, file="Rdata/decision_analytic/PSA_Only/mvglmer_psa.Rdata")

survModel = survreg(Surv(progression_time_start, progression_time_end, type = "interval2") ~ 
                      I(Age - 70) +  I((Age - 70)^2), data = prias.id, model = TRUE)
save(survModel, file="Rdata/decision_analytic/PSA_Only/survModel.Rdata")

forms = list("log2psaplus1" = "value",
             "log2psaplus1" = list(fixed = ~ 0 + dns(visitTimeYears, knots=c(0.1, 0.7, 4), Boundary.knots=c(0, 5.42)),
                                        random=~0 + dns(visitTimeYears, knots=c(0.1, 0.7, 4), Boundary.knots=c(0, 5.42)),
                                        indFixed = 4:7, indRandom=2:5, name = "slope"),
             "log2psaplus1" = list(fixed = ~ 0 + d2ns(visitTimeYears, knots=c(0.1, 0.7, 4), Boundary.knots=c(0, 5.42)),
                                   random=~0 + d2ns(visitTimeYears, knots=c(0.1, 0.7, 4), Boundary.knots=c(0, 5.42)),
                                   indFixed = 4:7, indRandom=2:5, name = "acceleration"))

startTime_mvJoint_psa = Sys.time()
mvJoint_psa = mvJointModelBayes(mvglmer_psa, survModel, timeVar = "visitTimeYears", Formulas = forms)
endTime_mvJoint_psa = Sys.time()
save(mvJoint_psa, file="Rdata/decision_analytic/PSA_Only/mvJoint_psa_accel.Rdata")
