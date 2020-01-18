library(JMbayes)
library(splines)
library(ggplot2)
library(ggpubr)

load("Rdata/gap3/PRIAS_2019/mvJoint_psa_time_scaled_light.Rdata")
load("Rdata/gap3/PRIAS_2019/cleandata.Rdata")
source("src/clinical_gap3/prediction_only_psa.R")

set.seed(2019)
pat_data = prias_long_final[prias_long_final$P_ID==102,]
##I am perturbing the PSA of one of the patients to demo effect of rising PSA
pat_data$log2psaplus1[c(nrow(pat_data)-1, nrow(pat_data))] =  pat_data$log2psaplus1[c(nrow(pat_data)-1, nrow(pat_data))] + runif(n = 2, 0.1, 0.5)
pat_data = pat_data[pat_data$year_visit<=3.0,]
pat_data$year_visit[nrow(pat_data)] = 2.5

#####################
# Step 1:
# First I will get the whole profile of the patient
# grid of future visits defined by time_lookahead_decision
#####################
time_lookahead_decision = seq(2.5, 10, 0.5)
M=500
latest_survival_time = 1.5
surv_lookahead_decision = getExpectedFutureOutcomes(object = mvJoint_psa_time_scaled,
                                                    patient_data = pat_data,
                                                    latest_survival_time = latest_survival_time,
                                                    earliest_failure_time = Inf,
                                                    survival_predict_times = time_lookahead_decision,
                                                    M = M)$predicted_surv_prob
mean_surv_lookahead_decision = rowMeans(surv_lookahead_decision)
#####################
# Step 2:
# Now we use a threshold on survival probability, and create a schedule,
# a minimum gap of 1 year b/w consecutive biopsies
#####################
surv_threshold = 0.9

proposed_biopsy_times = c()
latest_biopsy_time = latest_survival_time

surv_lookahead_decision_temp = surv_lookahead_decision
for(i in 1:length(time_lookahead_decision)){
  if(mean(surv_lookahead_decision_temp[i,], na.rm = T) <= surv_threshold){
    latest_biopsy_time = time_lookahead_decision[i]
    proposed_biopsy_times = c(proposed_biopsy_times, time_lookahead_decision[i])
    surv_lookahead_decision_temp = apply(surv_lookahead_decision, 2, FUN = function(x){x/x[i]})
  }
}

#and a compulsory biopsy at last lookahead time
proposed_biopsy_times = unique(c(proposed_biopsy_times, max(time_lookahead_decision)) )

print(proposed_biopsy_times)

#####################
# Step 3:
# Calculate expected number of biopsies using interval method
#####################
if(length(proposed_biopsy_times)==1){
  biop_intervals = list(c(latest_survival_time, proposed_biopsy_times))
}else{
  biop_intervals = c(latest_survival_time,
                     rep(proposed_biopsy_times[-length(proposed_biopsy_times)],each=2),
                     proposed_biopsy_times[length(proposed_biopsy_times)])
  biop_intervals = split(biop_intervals, rep(1:(length(biop_intervals)/2), each=2))
}

# then assume the mean survival probability to be the TRUTH (ignore MCMC)
mean_cumrisk_lookahead_decision = c(0, 1 - mean_surv_lookahead_decision)
names(mean_cumrisk_lookahead_decision)[1] = latest_survival_time

expected_number_biopsies = 0
for(j in 1:length(biop_intervals)){
  interval = as.character(biop_intervals[[j]])
  cum_risk_interval = diff(mean_cumrisk_lookahead_decision[interval])
  expected_number_biopsies = expected_number_biopsies + j * cum_risk_interval
}

expected_number_biopsies = expected_number_biopsies + j * mean_surv_lookahead_decision[as.character(biop_intervals[[j]][2])]
print(expected_number_biopsies)
#####################
# Step 4:
# Calculate expected number of biopsies using Dimitris' method
#####################
# assume the mean survival probability to be the TRUTH (ignore MCMC)
prob_of_surv_curve_being_truth = 1
expected_number_biopsies_dimitris = 0
for(i in 1:length(proposed_biopsy_times)){
  expected_number_biopsies_dimitris = expected_number_biopsies_dimitris + prob_of_surv_curve_being_truth
  if(i>1){
    prob_of_surv_curve_being_truth = mean_surv_lookahead_decision[as.character(proposed_biopsy_times[i])]
  }else{
    prob_of_surv_curve_being_truth = mean_surv_lookahead_decision[as.character(proposed_biopsy_times[i])]
  }
}
