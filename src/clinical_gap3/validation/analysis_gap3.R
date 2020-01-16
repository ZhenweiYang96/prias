cohort = "UCSF"

library(JMbayes)
library(splines)
library(survival)

longds = get(paste0("longdata_", cohort))
longds$age = longds$Age
longds$year_visit = longds$visitTimeYears
longds = droplevels(longds[!is.na(longds$log2psaplus1),])

longds.id = get(paste0("longdata_", cohort, ".id"))
longds.id$age = longds.id$Age
longds.id$latest_survival_time[longds.id$latest_survival_time==0] = 0.001

longds.id = longds.id[longds.id$P_ID %in% unique(longds$P_ID),]

survModel = survreg(Surv(latest_survival_time, earliest_failure_time, type = "interval2") ~ 
                      age, x = T, data = longds.id, model = TRUE, control = survreg.control(maxiter = 500))

startTime_mvglmer = Sys.time()
mvglmer_psa_time_scaled = mvglmer(list(log2psaplus1 ~ age +
                                         ns(I((year_visit-2)/2), knots=(c(0.5, 1.3, 3)-2)/2, Boundary.knots=(c(0, 6.3)-2)/2) + 
                                         (ns(I((year_visit-2)/2), knots=(c(0.5, 1.3, 3)-2)/2, Boundary.knots=(c(0, 6.3)-2)/2)|P_ID)),
                                  data=longds, families = list(gaussian), engine = "JAGS",
                                  control = list(n.iter=60000, n.burnin = 10000))
endTime_mvglmer = Sys.time()
print(endTime_mvglmer - startTime_mvglmer)

forms_psa = list("log2psaplus1" = "value",
                 "log2psaplus1" = list(fixed = ~ 0 + dns(I((year_visit-2)/2), knots=(c(0.5, 1.3, 3)-2)/2, Boundary.knots=(c(0, 6.3)-2)/2),
                                       random=~0 + dns(I((year_visit-2)/2), knots=(c(0.5, 1.3, 3)-2)/2, Boundary.knots=(c(0, 6.3)-2)/2),
                                       indFixed = 3:6, indRandom=2:5, name = "slope"))

startTime_mvJoint = Sys.time()
mvJoint_psa_time_scaled = mvJointModelBayes(mvglmer_psa_time_scaled, survModel, 
                                            timeVar = "year_visit", 
                                            Formulas = forms_psa, control = list(n_cores=3))
endTime_mvJoint = Sys.time()
