library(JMbayes)
library(splines)
library(ggplot2)
library(ggpubr)

source("src/clinical_gap3/prediction_only_psa.R")
load("Rdata/gap3/PRIAS_2019/mvJoint_psa_time_scaled_light.Rdata")
load("Rdata/gap3/PRIAS_2019/cleandata.Rdata")

set.seed(2019)
pat_data = prias_long_final[prias_long_final$P_ID==102,]
##I am perturbing the PSA of one of the patients to demo effect of rising PSA
pat_data$log2psaplus1[c(nrow(pat_data)-1, nrow(pat_data))] =  pat_data$log2psaplus1[c(nrow(pat_data)-1, nrow(pat_data))] + runif(n = 2, 0, 0.5)

pat_data = pat_data[pat_data$year_visit <= 4.5,]

M=750

set.seed(2019)
pred_times_1 = seq(2, 8, length.out = 200)
surv_prob_1 = getExpectedFutureOutcomes(mvJoint_psa_time_scaled,
                          pat_data, latest_survival_time = 2, 
                          survival_predict_times=pred_times_1[-1], 
                          psaDist = "Tdist", M = M)
surv_prob_1 = rbind(rep(1, M), surv_prob_1$predicted_surv_prob)


pred_times_2 = pred_times_1[pred_times_1>=4]
surv_prob_2 = getExpectedFutureOutcomes(mvJoint_psa_time_scaled,
                                        pat_data, latest_survival_time = 4, 
                                        survival_predict_times=pred_times_2[-1], 
                                        psaDist = "Tdist", M = M)
surv_prob_2 = rbind(rep(1, M), surv_prob_2$predicted_surv_prob)

pred_times_3 = pred_times_2
surv_prob_3 = surv_prob_1[pred_times_1 %in% pred_times_3,]
surv_prob_3 = apply(surv_prob_3, 2, FUN = function(x){x/x[1]})

ggplot() + 
  geom_line(aes(x=pred_times_1, y=rowMeans(surv_prob_1), color='T > 2')) +
  geom_line(aes(x=pred_times_2, y=rowMeans(surv_prob_2), color='T > 4, updated b')) +
  geom_line(aes(x=pred_times_3, y=rowMeans(surv_prob_3), color='T > 4, using b from T >2')) +
  theme_bw() +
  theme(legend.position = "bottom",
        text = element_text(size=12)) +
  scale_x_continuous(breaks = 2:8, limits = c(2,8)) +
  scale_y_continuous(breaks=seq(0,1, by = 0.25),
                     labels = paste(seq(0,1, by = 0.25)*100, "%"),
                     limits = c(0,1)) +
  xlab("Follow-up time (years)") + 
  ylab("Survival probability (%)")


ggplot() + 
  geom_line(aes(x=pred_times_2, y=rowMeans(surv_prob_2), color='updated b')) +
  geom_line(aes(x=pred_times_3, y=rowMeans(surv_prob_3), color='using b from T >2')) +
  geom_ribbon(aes(x=pred_times_2, ymin=apply(surv_prob_2,1, quantile, probs=0.025), 
                  ymax=apply(surv_prob_2,1, quantile, probs=0.975), fill='updated b'),alpha=0.15) +
  geom_ribbon(aes(x=pred_times_3, ymin=apply(surv_prob_3,1, quantile, probs=0.025), 
                  ymax=apply(surv_prob_3,1, quantile, probs=0.975), fill='using b from T >2'),alpha=0.15) +
  theme_bw() +
  theme(legend.position = "bottom",
        text = element_text(size=12)) +
  scale_x_continuous(breaks = 4:8, limits = c(4,8)) +
  scale_y_continuous(breaks=seq(0,1, by = 0.25),
                     labels = paste(seq(0,1, by = 0.25)*100, "%"),
                     limits = c(0,1)) +
  xlab("Follow-up time (years)") + 
  ylab("Survival probability (%)")
