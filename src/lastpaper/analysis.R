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
save(survModel, file="Rdata/gap3/PRIAS_2019/survModel.Rdata", version = 2)

prias_psa_dre = droplevels(prias_long_final[!(is.na(prias_long_final$psa) & is.na(prias_long_final$palpable_dre)),])
startTime_mvglmer = Sys.time()
mvglmer_psa_dre_time_scaled = mvglmer(list(palpable_dre ~ age + I((year_visit-2)/2)  + 
                                         (I((year_visit-2)/2)|P_ID),
                                       
                                       log2psaplus1 ~ age +
                                         ns(I((year_visit-2)/2), knots=(c(0.5, 1.3, 3)-2)/2, Boundary.knots=(c(0, 6.3)-2)/2) + 
                                         (ns(I((year_visit-2)/2), knots=(c(0.5, 1.3, 3)-2)/2, Boundary.knots=(c(0, 6.3)-2)/2)|P_ID)),
                                  data=prias_psa_dre, families = list(binomial, gaussian), engine = "STAN")

mvglmer_dre_psa = mvglmer(list(palpable_dre~age + year_visit  + 
                                 (year_visit|P_ID),
                               
                               log2psaplus1 ~ age + 
                                 ns(year_visit, knots=c(0.5, 1.3, 3), Boundary.knots=c(0, 6.3)) + 
                                 (ns(year_visit, knots=c(0.5, 1.3, 3), Boundary.knots=c(0, 6.3))|P_ID)),
                          
                          data=prias_psa_dre, families = list(binomial, gaussian), engine = "STAN",
                          control = list(n.iter=10000))

endTime_mvglmer = Sys.time()
print(endTime_mvglmer - startTime_mvglmer)

save(mvglmer_psa_dre_time_scaled, file="Rdata/gap3/PRIAS_2019/mvglmer_psa_dre_time_scaled.Rdata", version = 2)

forms_psa = list("log2psaplus1" = "value",
                 "log2psaplus1" = list(fixed = ~ 0 + dns(I((year_visit-2)/2), knots=(c(0.5, 1.3, 3)-2)/2, Boundary.knots=(c(0, 6.3)-2)/2),
                                       random=~0 + dns(I((year_visit-2)/2), knots=(c(0.5, 1.3, 3)-2)/2, Boundary.knots=(c(0, 6.3)-2)/2),
                                       indFixed = 3:6, indRandom=2:5, name = "slope"))

startTime_mvJoint = Sys.time()
mvJoint_psa_time_scaled = mvJointModelBayes(mvglmer_psa_time_scaled, survModel, 
                                            timeVar = "year_visit", 
                                            Formulas = forms_psa, control = list(n_cores=3))
endTime_mvJoint = Sys.time()

save(mvJoint_psa_time_scaled, file="Rdata/gap3/PRIAS_2019/mvJoint_psa_time_scaled.Rdata", version = 2)
