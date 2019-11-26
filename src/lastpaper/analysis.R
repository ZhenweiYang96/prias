load("Rdata/gap3/PRIAS_2019/cleandata.Rdata")

prias_final.id$latest_survival_time[prias_final.id$latest_survival_time==0] = 0.0001
prias_psa_dre = droplevels(prias_long_final[!(is.na(prias_long_final$psa) & is.na(prias_long_final$palpable_dre)),])

N_Cores = 6

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

##########################################
# Non parameteric model
#interval censoring NPMLE. Takes time to run. beware
##########################################
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


########################################
# Model 1:  2 internal knots, quadratic age
#######################################
startTime_mvglmer = Sys.time()
mvglmer_dre_psa_2knots_quad_age = mvglmer(list(palpable_dre~I(age-65) +  I((age-65)^2) + year_visit  +  (year_visit|P_ID),
                                               
                                               log2psaplus1~I(age-65) +  I((age-65)^2) + 
                                                 ns(year_visit, knots=c(0.75, 2.12), Boundary.knots=c(0, 6.4)) + 
                                                 (ns(year_visit, knots=c(0.75, 2.12), Boundary.knots=c(0, 6.4))|P_ID)),
                                          
                                          data=prias_psa_dre, families = list(binomial, gaussian), engine = "STAN",
                                          control = list(n.iter=10000, seed=2019))

endTime_mvglmer = Sys.time()
print(endTime_mvglmer - startTime_mvglmer)
save(mvglmer_dre_psa_2knots_quad_age, file="Rdata/lastpaper/fitted_model/mvglmer_2knots_quadage.Rdata")

survModel = survreg(Surv(latest_survival_time, earliest_failure_time, type = "interval2") ~ 
                      I(age-65) +  I((age-65)^2), x = T, data = prias_final.id, model = TRUE)

forms_dre_psa_2knots_quad_age = list("palpable_dre"= "value",
                                     "log2psaplus1" = "value",
                                     "log2psaplus1" = list(fixed = ~ 0 + dns(year_visit, knots=c(0.75, 2.12),Boundary.knots=c(0, 6.4)),
                                                           random=~0 + dns(year_visit, knots=c(0.75, 2.12), Boundary.knots=c(0, 6.4)),
                                                           indFixed = 4:6, indRandom=2:4, name = "slope"))

startTime_mvJoint = Sys.time()
mvJoint_dre_psa_2knots_quad_age = mvJointModelBayes(mvglmer_dre_psa_2knots_quad_age, survModel, 
                                                    timeVar = "year_visit", 
                                                    Formulas = forms_dre_psa_2knots_quad_age, 
                                                    control = list(n_cores=N_Cores, seed=2019))
endTime_mvJoint = Sys.time()
print(endTime_mvJoint - startTime_mvJoint)
save(mvJoint_dre_psa_2knots_quad_age, file="Rdata/lastpaper/fitted_model/mvJoint_dre_psa_2knots_quad_age.Rdata")


########################################
# Model 2:  2 internal knots, linear age
#######################################
startTime_mvglmer = Sys.time()
mvglmer_dre_psa_2knots_lin_age = mvglmer(list(palpable_dre~I(age-65) + year_visit  +  (year_visit|P_ID),
                                              
                                              log2psaplus1~I(age-65) + 
                                                ns(year_visit, knots=c(0.75, 2.12), Boundary.knots=c(0, 6.4)) + 
                                                (ns(year_visit, knots=c(0.75, 2.12), Boundary.knots=c(0, 6.4))|P_ID)),
                                         
                                         data=prias_psa_dre, families = list(binomial, gaussian), engine = "STAN",
                                         control = list(n.iter=10000, seed=2019))

endTime_mvglmer = Sys.time()
print(endTime_mvglmer - startTime_mvglmer)
save(mvglmer_dre_psa_2knots_lin_age, file="Rdata/lastpaper/fitted_model/mvglmer_2knots_lin_age.Rdata")

survModel = survreg(Surv(latest_survival_time, earliest_failure_time, type = "interval2") ~ I(age-65), x = T, data = prias_final.id, model = TRUE)

forms_dre_psa_2knots_lin_age = list("palpable_dre"= "value",
                                    "log2psaplus1" = "value",
                                    "log2psaplus1" = list(fixed = ~0 + dns(year_visit, knots=c(0.75, 2.12),Boundary.knots=c(0, 6.4)),
                                                          random = ~0 + dns(year_visit, knots=c(0.75, 2.12), Boundary.knots=c(0, 6.4)),
                                                          indFixed = 3:5, indRandom=2:4, name = "slope"))

startTime_mvJoint = Sys.time()
mvJoint_dre_psa_2knots_lin_age = mvJointModelBayes(mvglmer_dre_psa_2knots_lin_age, survModel, 
                                                   timeVar = "year_visit", 
                                                   Formulas = forms_dre_psa_2knots_lin_age, 
                                                   control = list(n_cores=N_Cores, seed=2019))
endTime_mvJoint = Sys.time()
print(endTime_mvJoint - startTime_mvJoint)
save(mvJoint_dre_psa_2knots_lin_age, file="Rdata/lastpaper/fitted_model/mvJoint_dre_psa_2knots_lin_age.Rdata")

########################################
# Model 3:  3 internal knots, quadratic age
#######################################
startTime_mvglmer = Sys.time()
mvglmer_dre_psa_3knots_quad_age = mvglmer(list(palpable_dre~I(age-65) +  I((age-65)^2) + year_visit  + 
                                                 (year_visit|P_ID),
                                               
                                               log2psaplus1 ~ I(age-65) +  I((age-65)^2) +
                                                 ns(year_visit, knots=c(0.5, 1.3, 3), Boundary.knots=c(0, 6.4)) + 
                                                 (ns(year_visit, knots=c(0.5, 1.3, 3), Boundary.knots=c(0, 6.4))|P_ID)),
                                          
                                          data=prias_psa_dre, families = list(binomial, gaussian), engine = "STAN",
                                          control = list(n.iter=10000, seed=2019))
endTime_mvglmer = Sys.time()
print(endTime_mvglmer - startTime_mvglmer)
save(mvglmer_dre_psa_3knots_quad_age, file="Rdata/lastpaper/fitted_model/mvglmer_3knots_quadage.Rdata")

survModel = survreg(Surv(latest_survival_time, earliest_failure_time, type = "interval2") ~ 
                      I(age-65) +  I((age-65)^2), x = T, data = prias_final.id, model = TRUE)

forms_dre_psa_3knots_quad_age = list("palpable_dre"= "value",
                                     "log2psaplus1" = "value",
                                     "log2psaplus1" = list(fixed = ~ 0 + dns(year_visit, knots=c(0.5, 1.3, 3),Boundary.knots=c(0, 6.4)),
                                                           random=~0 + dns(year_visit, knots=c(0.5, 1.3, 3), Boundary.knots=c(0, 6.4)),
                                                           indFixed = 4:7, indRandom=2:5, name = "slope"))

startTime_mvJoint = Sys.time()
mvJoint_dre_psa_3knots_quad_age= mvJointModelBayes(mvglmer_dre_psa_3knots_quad_age, survModel, 
                                                   timeVar = "year_visit", 
                                                   Formulas = forms_dre_psa_3knots_quad_age, 
                                                   control = list(n_cores=N_Cores, seed=2019))
endTime_mvJoint = Sys.time()
print(endTime_mvJoint - startTime_mvJoint)
save(mvJoint_dre_psa_3knots_quad_age, file="Rdata/lastpaper/fitted_model/mvJoint_dre_psa_3knots_quad_age.Rdata")