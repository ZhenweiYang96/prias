library(lme4)
library(splines)

source("common.R")

dre_data_set =  prias_long[!is.na(prias_long$didre),]
dre_data_set = prias_long_survival_eligible[!is.na(prias_long_survival_eligible$didre),]

dre_data_set$P_ID = droplevels(dre_data_set$P_ID)
dre_data_set$didre_num = ifelse(dre_data_set$didre=="T1", 0, 1)
dre_data_set$segmentVisitTimeYears =  (dre_data_set$visitTimeYears-1)*(dre_data_set$visitTimeYears>1)
save.image(file = "dre_gl.Rdata")

quantiles_time_years = round(quantile(dre_data_set$visitTimeYears, probs = seq(0,1, by=0.005)), 3)
quantiles_time_years = unique(quantiles_time_years)

odds = vector("numeric", length(quantiles_time_years))
quantiles_time_years = c(-1, quantiles_time_years)
for(i in 2:length(quantiles_time_years)){
  prias_df_temp = dre_data_set[dre_data_set$visitTimeYears<=quantiles_time_years[i] & 
                                 dre_data_set$visitTimeYears > quantiles_time_years[i-1],]
  print(paste(quantiles_time_years[i-1], quantiles_time_years[i], length(unique(prias_df_temp$P_ID))))
  odds[i-1] = table(prias_df_temp$didre)[">T1"]/table(prias_df_temp$didre)["T1"]
}

qplot(y=log(odds), x=quantiles_time_years[-1], geom=c("point","line", "smooth")) + 
  ticksX(from = 0, max = 10, by =  1) + xlim(0,7) + xlab("Time (years)") + ylab("log (Odds: >T1 vs T1))")

ggplot(data = dre_data_set, aes(y=didre_num,x=visitTimeYears)) + 
  geom_point() + stat_smooth(geom="line", colour="blue") + ticksX(from = 0, max = 10, 1) + xlim(0,10) + 
  xlab("Time (years)") + ylab("Probability (DRE > T1c)")


#Model fitting
mvglm_dre_segmented_jags = mvglmer(list(didre_num~Age + visitTimeYears +
                                              I((visitTimeYears-1)*(visitTimeYears>1)) + 
                                              (1|P_ID)),
                                       data = dre_data_set, families = list(binomial), engine="JAGS")
save.image(file = "dre_gl.Rdata")
mvJoint_dre_segmented_jags = mvJointModelBayes(mvglm_dre_segmented_jags, coxModel, timeVar = "visitTimeYears")
save.image(file = "dre_gl.Rdata")

###############################
#Some graphs based on these estimates
##############################
simVisitTime = seq(0, 10, 0.5)
longBetas = summary(mvJoint_dre_segmented_jags)$Outcome1$PostMean
varIntercept = summary(mvJoint_dre_segmented_jags)$D
logOddsMedianSub = longBetas[1] + longBetas[2] * 70 + longBetas[3] * simVisitTime + 
  longBetas[4] * (simVisitTime - 1) * (simVisitTime > 1)
probabilityMedianSub= exp(logOddsMedianSub) / (1+exp(logOddsMedianSub))
logOddsAvgSubject = logOddsMedianSub/sqrt(1+0.346 * varIntercept)
probabilityAvgSubject = exp(logOddsAvgSubject) / (1+exp(logOddsAvgSubject))

ggplot(data=data.frame(time=rep(simVisitTime,2), 
                       probability =c(probabilityMedianSub, probabilityAvgSubject), 
                       Type=c(rep("Median Subject",length(simVisitTime)), rep("Average Subject",length(simVisitTime))))) + 
  geom_line(aes(x=time, y=probability, group=Type, color=Type)) +
  ticksX(from = 0, max = 10, by = 0.5) + 
  ticksY(from = 0, max = 0.15, by = 0.01) + xlab("Follow-up time (Years)") + 
  ylab("Probability(DRE > T1c)") +  theme(legend.position = c(.85, .85))

sim_dre_prob_grt_T1c = seq(0.01, 0.01 * 10, by=0.005)
sim_dre_odds_grt_T1c = sim_dre_prob_grt_T1c / (1-sim_dre_prob_grt_T1c)
alpha = summary(mvJoint_dre_segmented_jags)$Survival$PostMean[3]
simhr = exp(alpha * (sim_dre_odds_grt_T1c - sim_dre_odds_grt_T1c[1]))

qplot(x = sim_dre_prob_grt_T1c, y=simhr, geom = "line", ylab = "Hazard Ratio", xlab = "Probability(DRE > T1c)") + 
  ticksX(from=sim_dre_prob_grt_T1c[1], max = max(sim_dre_prob_grt_T1c), by = 0.01) + 
  ticksY(from=0, 1.01, by = 0.001)

