load("Rdata/gap3/PRIAS_2019/cleandata.Rdata")

library(JMbayes)

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("Icens")
# install.packages("interval")

library(Icens)
library(interval)
library(splines)
library(survival)

#interval censoring NPMLE. Takes time to run. beware
npmle_prias = icfit(Surv(latest_survival_time,earliest_failure_time,type="interval2")~1, 
                    data=prias_final.id, control=icfitControl(maxit = 20000))
save(npmle_prias, file="Rdata/gap3/PRIAS_2019/npmle_prias.Rdata")
#The density function is given by npmle_prias$pf
survProb = 1 - cumsum(npmle_prias$pf)
survProb = c(1, survProb)

survIntervals = npmle_prias$intmap
survIntervals = cbind(c(0,0), survIntervals)

timePoints = as.numeric(survIntervals)
survProbs = c(1,as.numeric(rep(survProb, each=2)))[1:length(timePoints)]

FONT_SIZE=11
npmle_plot = ggplot() + geom_line(aes(x=timePoints, y=1-survProbs)) +  
  coord_cartesian(xlim=c(0,10)) + 
  theme_bw() +
  theme(text = element_text(size=FONT_SIZE), axis.text=element_text(size=FONT_SIZE),
        axis.line = element_line(), 
        plot.margin = margin(0, 0, 0, 0, "pt")) + 
  scale_y_continuous(breaks = seq(0, 1, 0.25), labels = paste0(seq(0, 1, 0.25)*100, "%"),
                     limits = c(0,1)) + 
  ylab("Cumulative risk of cancer progression (%)") +
  xlab("Follow-up time (years)")

## Now with joint model
prias_final.id$latest_survival_time[prias_final.id$latest_survival_time==0] = 0.0001
survModel = survreg(Surv(latest_survival_time, earliest_failure_time, type = "interval2") ~ 
                      age, x = T,
                    data = prias_final.id, model = TRUE)
save(survModel, file="Rdata/gap3/PRIAS_2019/survModel.Rdata")

model1 = lme(fixed=log2psaplus1 ~ age + ns(year_visit, knots=c(0.5, 1.3, 3), Boundary.knots=c(0, 6.3)), 
             random = ~ns(year_visit, knots=c(0.5, 1.3, 3), Boundary.knots=c(0, 6.3))|P_ID, 
             data = prias_psa, control = lmeControl(opt = "optim", optimMethod = "L-BFGS-B"))

prias_psa = droplevels(prias_long_final[!is.na(prias_long_final$psa),])
startTime_mvglmer = Sys.time()
mvglmer_psa = mvglmer(list(log2psaplus1 ~ age +
                                 ns(year_visit, knots=c(0.5, 1.3, 3), Boundary.knots=c(0, 6.3)) + 
                                 (ns(year_visit, knots=c(0.5, 1.3, 3), Boundary.knots=c(0, 6.3))|P_ID)),
                          data=prias_psa, families = list(gaussian), engine = "STAN",
                          control = list(n.iter=10000))
endTime_mvglmer = Sys.time()

print(endTime_mvglmer - startTime_mvglmer)
save(mvglmer_psa, file="Rdata/gap3/PRIAS_2019/mvglmer_psa_stanhome.Rdata")
save(mvglmer_psa, file="Rdata/gap3/PRIAS_2019/mvglmer_psa_stanpc2.Rdata")

forms_psa = list("log2psaplus1" = "value",
                 "log2psaplus1" = list(fixed = ~ 0 + dns(year_visit, knots=c(0.5, 1.3, 3), Boundary.knots=c(0, 6.3)),
                                       random=~0 + dns(year_visit, knots=c(0.5, 1.3, 3), Boundary.knots=c(0, 6.3)),
                                       indFixed = 3:6, indRandom=2:5, name = "slope"))

startTime_mvJoint = Sys.time()
mvJoint_psa = mvJointModelBayes(mvglmer_psa, survModel, 
                                timeVar = "year_visit", 
                                Formulas = forms_psa, control = list(n_cores=5))
endTime_mvJoint = Sys.time()

save(mvJoint_psa, file="Rdata/gap3/PRIAS_2019/mvJoint_psa_stanhome.Rdata")
save(mvJoint_psa, file="Rdata/gap3/PRIAS_2019/mvJoint_psa_stanpc2.Rdata")

rm(mvglmer_psa)
rm(mvJoint_psa)

################
# DRE PSA both
################
prias_psa_dre = droplevels(prias_long_final[!(is.na(prias_long_final$psa) & 
                                                is.na(prias_long_final$palpable_dre)),])
startTime_mvglmer = Sys.time()
mvglmer_dre_psa = mvglmer(list(palpable_dre~age + year_visit  + 
                                 (year_visit|P_ID),
                               
                               log2psaplus1 ~ age +
                                 ns(year_visit, knots=c(0.5, 1.3, 3), Boundary.knots=c(0, 6.3)) + 
                                 (ns(year_visit, knots=c(0.5, 1.3, 3), Boundary.knots=c(0, 6.3))|P_ID)),
                          data=prias_psa_dre, families = list(binomial, gaussian), 
                          engine = "STAN",
                          control = list(n.iter=10000))
endTime_mvglmer = Sys.time()
print(endTime_mvglmer - startTime_mvglmer)

save(mvglmer_dre_psa, file="Rdata/gap3/PRIAS_2019/mvglmer_dre_psa.Rdata")

forms_dre_value = list("palpable_dre" = "value",
                       "log2psaplus1" = "value",
                       "log2psaplus1" = list(fixed = ~ 0 + dns(year_visit, knots=c(0.5, 1.3, 3), Boundary.knots=c(0, 6.3)),
                                             random=~0 + dns(year_visit, knots=c(0.5, 1.3, 3), Boundary.knots=c(0, 6.3)),
                                             indFixed = 3:5, indRandom=2:5, name = "slope"))

startTime_mvJoint_dre_value = Sys.time()
mvJoint_dre_psa = mvJointModelBayes(mvglmer_dre_psa, survModel, 
                                              timeVar = "year_visit", 
                                              Formulas = forms_dre_value,
                                              control = list(n_cores=5))
endTime_mvJoint_dre_value = Sys.time()
print(endTime_mvJoint_dre_value - startTime_mvJoint_dre_value)
save(mvJoint_dre_psa, file="Rdata/gap3/PRIAS_2019/mvJoint_dre_psa.Rdata")
rm(mvJoint_dre_psa)

