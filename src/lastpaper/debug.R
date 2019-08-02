load("Rdata/lastpaper/fpr/auc_cumrisk.Rdata")
load("Rdata/lastpaper/fpr/aucrisk.Rdata")
fpr = combined[,-c(1:12)]
colnames(fpr) = paste0("fpr_", colnames(fpr))
rm(combined)

load("Rdata/lastpaper/fpr/threshold.Rdata")
threshold = testDs.id[,-c(1:12)]
colnames(threshold) = paste0("threshold_", colnames(threshold))
rm(testDs.id)

load("Rdata/lastpaper/sim_seed_2019_sched_wtbyhorizonriskT.Rdata")
delay_exchange_T = testDs.id[,-c(1:12)]
colnames(delay_exchange_T) = paste0("delayexchangeT_", colnames(delay_exchange_T))
rm(testDs.id)

load("Rdata/lastpaper/sim_seed_2019_sched_wtbyhorizonriskFALSE.Rdata")
delay_exchange_F = testDs.id[,-c(1:12)]
colnames(delay_exchange_F) = paste0("delayexchangeF_", colnames(delay_exchange_F))

testDs.id = cbind(testDs.id[,1:12], fpr, threshold, delay_exchange_T, delay_exchange_F)

rm(fpr, threshold, delay_exchange_T, delay_exchange_F)

methodNames = colnames(testDs.id)[seq(13, ncol(testDs.id), 2)]
testDs.id_long = reshape(testDs.id, varying = list(seq(13, ncol(testDs.id), 2), 
                                                   seq(14, ncol(testDs.id), 2)), 
                         v.names = c("nb", "delay"), idvar = "P_ID", direction = "long", 
                         timevar = "methodName",
                         times = methodNames)
testDs.id_long$methodtype = factor(c(rep("fpr", 21*300), rep("threshold", 21*300),
                                     rep("delay_exchangeT", 9 *300), rep("delay_exchangeF", 9 *300)))

testDs.id_long$methodName = factor(testDs.id_long$methodName)

temp_filter = is.na(testDs.id_long$delay) | testDs.id_long$delay < 0
testDs.id_long$nb[temp_filter] = testDs.id_long$nb[temp_filter] + 1
testDs.id_long$delay[temp_filter] = 10 - testDs.id_long$progression_time[temp_filter]

combined_summary = data.frame(methodName=levels(testDs.id_long$methodName),
                              median_nb = as.numeric(by(testDs.id_long$methodName, data = testDs.id_long$nb, median)),
                              q1_nb = as.numeric(by(testDs.id_long$methodName, data = testDs.id_long$nb, quantile, probs=0.25)),
                              q3_nb = as.numeric(by(testDs.id_long$methodName, data = testDs.id_long$nb, quantile, probs=0.75)),
                              mean_nb = as.numeric(by(testDs.id_long$methodName, data = testDs.id_long$nb, mean)),
                              median_nb_nonprog = as.numeric(by(testDs.id_long$methodName[testDs.id_long$progression_time==10], data = testDs.id_long$nb[testDs.id_long$progression_time==10], median)),
                              q1_nb_nonprog = as.numeric(by(testDs.id_long$methodName[testDs.id_long$progression_time==10], data = testDs.id_long$nb[testDs.id_long$progression_time==10], quantile, probs=0.25)),
                              q3_nb_nonprog = as.numeric(by(testDs.id_long$methodName[testDs.id_long$progression_time==10], data = testDs.id_long$nb[testDs.id_long$progression_time==10], quantile, probs=0.75)),
                              mean_nb_nonprog = as.numeric(by(testDs.id_long$methodName[testDs.id_long$progression_time==10], data = testDs.id_long$nb[testDs.id_long$progression_time==10], mean)),
                              median_nb_prog = as.numeric(by(testDs.id_long$methodName[testDs.id_long$progression_time<10], data = testDs.id_long$nb[testDs.id_long$progression_time<10], median)),
                              q1_nb_prog = as.numeric(by(testDs.id_long$methodName[testDs.id_long$progression_time<10], data = testDs.id_long$nb[testDs.id_long$progression_time<10], quantile, probs=0.25)),
                              q3_nb_prog = as.numeric(by(testDs.id_long$methodName[testDs.id_long$progression_time<10], data = testDs.id_long$nb[testDs.id_long$progression_time<10], quantile, probs=0.75)),
                              mean_nb_prog = as.numeric(by(testDs.id_long$methodName[testDs.id_long$progression_time<10], data = testDs.id_long$nb[testDs.id_long$progression_time<10], mean)),
                              median_delay_prog = as.numeric(by(testDs.id_long$methodName[testDs.id_long$progression_time<10], data = testDs.id_long$delay[testDs.id_long$progression_time<10], median)),
                              q1_delay_prog = as.numeric(by(testDs.id_long$methodName[testDs.id_long$progression_time<10], data = testDs.id_long$delay[testDs.id_long$progression_time<10], quantile, probs=0.25)),
                              q3_delay_prog = as.numeric(by(testDs.id_long$methodName[testDs.id_long$progression_time<10], data = testDs.id_long$delay[testDs.id_long$progression_time<10], quantile, probs=0.75)),
                              mean_delay_prog = as.numeric(by(testDs.id_long$methodName[testDs.id_long$progression_time<10], data = testDs.id_long$delay[testDs.id_long$progression_time<10], mean)))
combined_summary$methodtype = rep(levels(testDs.id_long$methodtype), table(testDs.id_long$methodtype)/nrow(testDs.id))

#nb all patients, but delay for progressions
ggplot(combined_summary) + 
  geom_point(aes(x=mean_nb, y=mean_delay_prog, color=methodtype)) + 
  theme_bw() + theme(legend.position = "bottom") + 
  scale_x_continuous(breaks=seq(1,10,1)) +
  scale_y_continuous(breaks=seq(1,7,1)) +
  facet_wrap(~methodtype) + 
  xlab("Number of biopsies") + ylab("Delay (years)")

#nb and delay for progressions
ggplot(combined_summary) + 
  geom_point(aes(x=median_nb_prog, y=median_delay_prog, color=methodtype)) + 
  theme_bw() + theme(legend.position = "bottom") + 
  scale_x_continuous(breaks=seq(1,10,1)) +
  scale_y_continuous(breaks=seq(1,7,1)) +
  xlab("Number of biopsies") + ylab("Delay (years)")

thres_filter = (combined_long$threshold>0.4 & combined_long$threshold<0.7) | combined_long$threshold>1
#nb boxplot non progressions
ggplot(combined_long[combined_long$progression_time==10 & thres_filter,]) + 
  geom_boxplot(aes(x=factor(threshold), y=nb, color=methodtype)) + 
  theme_bw() + theme(legend.position = "bottom") 

#nb boxplot progressions
ggplot(combined_long[combined_long$progression_time<10& thres_filter,]) + 
  geom_boxplot(aes(x=factor(threshold), y=nb, color=methodtype)) + 
  theme_bw() + theme(legend.position = "bottom") 

#delay boxplot progressions
ggplot(combined_long[combined_long$progression_time<10& thres_filter,]) + 
  geom_boxplot(aes(x=factor(threshold), y=delay, color=methodtype)) + 
  theme_bw() + theme(legend.position = "bottom") 
