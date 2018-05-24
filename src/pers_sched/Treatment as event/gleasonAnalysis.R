library(lme4)
library(splines)

source("common.R")

quantiles_time_years = round(quantile(prias_long[!is.na(prias_long$gleason),]$visitTimeYears, probs = seq(0,1, by=0.005)), 3)
quantiles_time_years = unique(quantiles_time_years)

odds = vector("numeric", length(quantiles_time_years))
quantiles_time_years = c(-1, quantiles_time_years)
for(i in 2:length(quantiles_time_years)){
  prias_df_temp = prias_long[prias_long$visitTimeYears<=quantiles_time_years[i] & 
                               prias_long$visitTimeYears > quantiles_time_years[i-1],]
  print(paste(quantiles_time_years[i-1], quantiles_time_years[i], length(unique(prias_df_temp$P_ID))))
  odds[i-1] = table(prias_df_temp$digleason)["High"]/table(prias_df_temp$digleason)["Low"]
}

qplot(y=log(odds), x=quantiles_time_years[-1], geom=c("point","line", "smooth")) + 
  ticksX(0, 15, 1) + xlab("Time (years)") + ylab("log(Odds: High vs. Low)")

ggplot(data = prias_long, aes(y=(as.numeric(digleason)-1),x=visitTimeYears)) + 
  geom_point() + stat_smooth(geom="line", colour="blue") + ticksX(from = 0, max = 10, 1) + xlim(0,10) + 
  xlab("Time (years)") + ylab("Probability (Gleason = High)")

#Model fitting
gleason_data_set =  prias_long[!is.na(prias_long$digleason),]
gleason_data_set$P_ID = droplevels(gleason_data_set$P_ID)

#There are ~8000 obs and ~5000 ppl. So I am using only random intercept. however in the fixed effects
#I ti

gleason_data_set$firstvisit = ifelse(gleason_data_set$visit_number==1, 1, 0)
gleason_data_set$digleason_num = ifelse(gleason_data_set$digleason=="Low", 0, 1)
gleason_data_set$segmentVisitTimeYears =  (gleason_data_set$visitTimeYears-1)*(gleason_data_set$visitTimeYears>1)

write.csv(gleason_data_set, file = "gleason_long.csv", row.names = F)

#Doesn't work
# glmer_gleason_int = glmer(digleason ~I(Age/10) + I(10 * firstvisit) + visitTimeYears + 
#                             (1|P_ID), data = gleason_data_set, family = binomial)
# 
# glmer_gleason_int = glmer(digleason ~I(Age/10) + ns(visitTimeYears, knots=c(0, 0.1), Boundary.knots=c(0,1.05)) +
#         (1|P_ID), data = gleason_data_set, family = binomial)
# 
# 
glmer_gleason_int_glmmpql = glmmPQL(digleason_num ~  Age + segmentVisitTimeYears+ visitTimeYears, random = ~ 1 | P_ID,
         family = "binomial", data = gleason_data_set)
 
dLongBin <- function (y, eta.y, scale, log = FALSE, data) {
  dbinom(x = y, size = 1, prob = plogis(eta.y), log = log)
}

jointFit_gleason <- jointModelBayes(glmer_gleason_int_glmmpql, coxModel, timeVar = "visitTimeYears",
                                 densLong = dLongBin)

mvglm_gleason_stan = mvglmer(list(digleason_num~Age + visitTimeYears + (1|P_ID)),
                        data = gleason_data_set, families = list(binomial), engine="STAN")

mvglm_gleason_jags = mvglmer(list(digleason_num~Age + visitTimeYears + (1|P_ID)),
                        data = gleason_data_set, families = list(binomial), engine="JAGS")

save.image(file = "gleason.Rdata")

mvglm_gleason_firstvisit_jags = mvglmer(list(digleason_num~Age + firstvisit + visitTimeYears + (1|P_ID)),
                             data = gleason_data_set, families = list(binomial), engine="JAGS")


mvglm_gleason_segmented_jags = mvglmer(list(digleason_num~Age + visitTimeYears +
                                              I((visitTimeYears-0.05)*(visitTimeYears>0.05)) + 
                                             (1|P_ID)),
                                        data = gleason_data_set, families = list(binomial), engine="JAGS")



jointmodel_gleason_jags = mvJointModelBayes(mvglm_gleason_jags, coxModel, timeVar = "visitTimeYears")
jointmodel_gleason_firstvisit_jags = mvJointModelBayes(mvglm_gleason_firstvisit_jags, coxModel, timeVar = "visitTimeYears")

jointmodel_gleason_segmented_jags = mvJointModelBayes(mvglm_gleason_segmented_jags, coxModel, timeVar = "visitTimeYears")
                  # Formulas = list("digleason_num" = "value",
                  #                 "digleason_num" = list(fixed = ~ 0 + dns(visitTimeYears, knots=c(0.1), Boundary.knots=c(0,1.05)),
                  #                                    random=~0,
                  #                                    indFixed = 3:4, name = "slope")))



###############################
#Some graphs based on these estimates
##############################
simVisitTime = seq(0,10, 0.5)
logOddsMedianSub = -7.9815 + 0.0446 * 70 + 6.40853 * simVisitTime - 6.33368 * (simVisitTime - 0.5) * (simVisitTime > 0.5)
probabilityMedianSub= exp(logOddsMedianSub) / (1+exp(logOddsMedianSub))
logOddsAvgSubject = logOddsMedianSub/sqrt(1+0.346*0.829*0.829)
probabilityAvgSubject = exp(logOddsAvgSubject) / (1+exp(logOddsAvgSubject))

ggplot(data=data.frame(time=rep(simVisitTime,2), 
                       probability =c(probabilityMedianSub, probabilityAvgSubject), 
                       Type=c(rep("Median Subject",length(simVisitTime)), rep("Average Subject",length(simVisitTime))))) + 
  geom_line(aes(x=time, y=probability, group=Type, color=Type)) +
  ticksX(from = 0, max = 10, by = 0.5) + 
 xlab("Follow-up time (Years)") + 
  ylab("Probability(Gleason = High)") +  theme(legend.position = c(.85, .15))

sim_gleason_prob_high = seq(0.01, 0.01 * 30, by=0.005)
sim_gleason_odds_high = sim_gleason_prob_high / (1-sim_gleason_prob_high)
simhr = exp(2.1131 * (sim_gleason_odds_high - sim_gleason_odds_high[1]))

qplot(x = sim_gleason_prob_high, y=simhr, geom = "line", ylab = "Hazard Ratio", xlab = "Probability(Gleason = High)") + 
  ticksX(from=sim_gleason_prob_high[1], max = max(sim_gleason_prob_high), by = 0.02) + 
  ticksY(from=1, max=3.5, by = 0.1)
