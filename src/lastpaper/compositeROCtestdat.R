newtestdata = sim_res$testData$testDs[sim_res$testData$testDs$progression_time>1,]

set.seed(2019)
cum_risk_T_start_onwards_test = by(INDICES = newtestdata$P_ID, data = newtestdata, FUN = function(pat_data){
  rowMeans(1 - getExpectedFutureOutcomes(fitted_jm, pat_data, 
                                         latest_survival_time = 1,
                                         earliest_failure_time = Inf,
                                         survival_predict_times = BIOPSY_TEST_TIMES[BIOPSY_TEST_TIMES >= 2],
                                         psaDist = "Tdist")$predicted_surv_prob)
})

decision = sapply(cum_risk_T_start_onwards_test, "[",1) > 0.1
real = newtestdata$progression_time[!duplicated(newtestdata$P_ID)] <=2

table(decision, real)

#     real
# decision FALSE TRUE
# FALSE   141   12
# TRUE     91   23

# versus
decision = sapply(cum_risk_T_start_onwards_test, "[",1) > 0.092 & 
  sapply(cum_risk_T_start_onwards_test, "[",3) > 0.168
real = newtestdata$progression_time[!duplicated(newtestdata$P_ID)] <=2

table(decision, real)
#     real
# decision FALSE TRUE
# FALSE   147   10
# TRUE     85   25

# versus
decision = sapply(cum_risk_T_start_onwards_test, "[",1) > 0.089 & 
  sapply(cum_risk_T_start_onwards_test, "[",3) > 0.169 &
  sapply(cum_risk_T_start_onwards_test, "[",5) > 0.235
real = newtestdata$progression_time[!duplicated(newtestdata$P_ID)] <=2

table(decision, real)
