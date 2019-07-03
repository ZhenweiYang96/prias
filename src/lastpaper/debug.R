T_start=2
T_horiz=3
horizon=10
M=0

#First make sure that longitudinal data is available only upto T_horiz
newdata = prias_long_final[prias_long_final$year_visit <= T_horiz,]

#Now select patients who have a progression time more than T_start
#This automatically means that they had no event until T_start
newdata = newdata[newdata$right_cens_time > T_start,]
newdata.id = newdata[!duplicated(newdata$P_ID),]

########################
#Make sure there are patients to calculate TPR etc
########################


#Now lets start with real patient status
newdata.id$real_period_status = rep(NA, nrow(newdata.id))
newdata.id$real_period_status[newdata.id$reclassification==T & newdata.id$right_cens_time <= T_horiz] = 1
newdata.id$real_period_status[newdata.id$reclassification==T & newdata.id$right_cens_time > T_horiz] = 0
newdata.id$real_period_status[newdata.id$reclassification==F & newdata.id$right_cens_time >= T_horiz] = 0

#Now for patients who had reclassification = F and right_cens_time <= T_horiz
# Their latest biopsy time is same right_cens_time
# we need to calculate their cum risk using their real data
subset_patients = newdata[newdata$reclassification==F & newdata$right_cens_time < T_horiz,]
subset_patients_cum_risk_Thoriz = by(INDICES = subset_patients$P_ID, data = subset_patients, FUN = function(pat_data){
  rowMeans(1 - getExpectedFutureOutcomes(mvJoint_psa_time_scaled, pat_data, 
                                         latest_survival_time = pat_data$right_cens_time[1],
                                         earliest_failure_time = Inf,
                                         survival_predict_times = T_horiz,
                                         psaDist = "Tdist")$predicted_surv_prob)
})

newdata.id$real_period_status[newdata.id$reclassification==F & 
                                newdata.id$right_cens_time < T_horiz] = as.numeric(subset_patients_cum_risk_Thoriz)

#Decision Epochs in years
PSA_CHECK_UP_TIME = c(seq(0, 2, 0.25), seq(2.5, horizon, 0.5))
#DRE check up time years
DRE_CHECK_UP_TIME = seq(0, horizon, 0.5)
BIOPSY_TEST_TIMES = DRE_CHECK_UP_TIME

#Now for each of the patients obtain their cumulative risk at 
#various follow-ups under the condition
#that they didnt have an event until T_start
cum_risk_T_start_onwards = by(INDICES = newdata$P_ID, data = newdata, FUN = function(pat_data){
  rowMeans(1 - getExpectedFutureOutcomes(mvJoint_psa_time_scaled, pat_data, 
                                         latest_survival_time = T_start,
                                         earliest_failure_time = Inf,
                                         survival_predict_times = BIOPSY_TEST_TIMES[BIOPSY_TEST_TIMES > T_start],
                                         psaDist = "Tdist")$predicted_surv_prob)
})


########Simple rule results
cum_risk_T_start_T_horiz = sapply(cum_risk_T_start_onwards, "[", as.character(T_horiz))
thresholds = c(seq(0,0.3, length.out = 1e6), seq(0.3, 1, length.out = 100))
results = t(sapply(thresholds, FUN = function(threshold){
  predicted_cancer = (cum_risk_T_start_T_horiz >= threshold)
  real_cancer = newdata.id$real_period_status
  
  c("threshold" = threshold,
    "nTP" = sum(real_cancer * predicted_cancer),
    "nFN" = sum(real_cancer * (1-predicted_cancer)),
    "nFP" = sum((1-real_cancer) * predicted_cancer),
    "nTN" = sum((1-real_cancer) * (1-predicted_cancer)))
}))

results = data.frame(results)
results$tpr = results$nTP / (results$nTP + results$nFN)
results$fpr = results$nFP / (results$nTN + results$nFP)


ggplot() + 
  geom_line(aes(x=results$fpr, y=results$tpr, color='simple')) + 
  geom_abline(intercept = 0, slope=1) +
  xlim(0,1) + ylim(0,1) + theme(legend.position = "bottom")


########## Composite rule results
cum_risk_T_start_T_horiz = sapply(cum_risk_T_start_onwards, "[", as.character(T_horiz))
cum_risk_T_start_T_horizplus1 = sapply(cum_risk_T_start_onwards, "[", as.character(T_horiz+1))
cum_risk_T_start_T_horizplus2 = sapply(cum_risk_T_start_onwards, "[", as.character(T_horiz+2))
cum_risk_T_start_T_horizplus3 = sapply(cum_risk_T_start_onwards, "[", as.character(T_horiz+3))

threshold_comp = expand.grid(thresholds[thresholds<=0.1], 
                             thresholds[thresholds<=0.2], 
                             thresholds[thresholds<=0.3], 
                             thresholds[thresholds<=0.4])

results_comp = t(sapply(1:nrow(threshold_comp), FUN = function(i){
  threshold1 = threshold_comp[i,1]
  threshold2 = threshold_comp[i,2]
  threshold3 = threshold_comp[i,3]
  threshold4 = threshold_comp[i,4]
  
  predicted_cancer = (cum_risk_T_start_T_horiz >= threshold1) &
    (cum_risk_T_start_T_horizplus1 >= threshold2) &
    (cum_risk_T_start_T_horizplus2 >= threshold3) &
    (cum_risk_T_start_T_horizplus3 >= threshold4)
  
  real_cancer = newdata.id$real_period_status
  
  c("threshold1" = threshold1,
    "threshold2" = threshold2,
    "threshold3" = threshold3,
    "threshold4" = threshold4,
    "nTP" = sum(real_cancer * predicted_cancer),
    "nFN" = sum(real_cancer * (1-predicted_cancer)),
    "nFP" = sum((1-real_cancer) * predicted_cancer),
    "nTN" = sum((1-real_cancer) * (1-predicted_cancer)))
}))

results_comp = data.frame(results_comp)
results_comp$tpr = results_comp$nTP / (results_comp$nTP + results_comp$nFN)
results_comp$fpr = results_comp$nFP / (results_comp$nTN + results_comp$nFP)

ggplot() + geom_line(aes(x=results_comp$fpr, y=results_comp$tpr, color='comp')) + 
  geom_line(aes(x=results$fpr, y=results$tpr, color='simple')) + 
  geom_abline(intercept = 0, slope=1) +
  xlim(0,1) + ylim(0,1) + theme(legend.position = "bottom")
