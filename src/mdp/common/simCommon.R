#Select parameters common for all types of algorithms
NR_DISCRETIZED_OBS = 100
MIN_BIOPSY_GAP = 1

N_MCMC_ITER = 200

#Constants
MAX_FOLLOW_UP_TIME = 10
#Decision Epochs in years
PSA_CHECK_UP_TIME = c(seq(0, 2, 0.25), seq(2.5, MAX_FOLLOW_UP_TIME, 0.5))
#DRE check up time years
DRE_CHECK_UP_TIME = c(seq(0, MAX_FOLLOW_UP_TIME, 0.5))

#Actions
BIOPSY = "B"
WAIT = "W"
ACTION_VECTOR = c(BIOPSY, WAIT)

#States
G6 = "G6"
G7 = "G7"
AT = "AT"
STATE_VECTOR = c(G6, G7, AT)

DISCOUNT_FACTOR = 1
DISCOUNT_FACTORS = DISCOUNT_FACTOR^(1:length(PSA_CHECK_UP_TIME))
names(DISCOUNT_FACTORS) = PSA_CHECK_UP_TIME

source("src/mdp/common/prediction.R")

#########################################
# We define various functions from here onwards
#########################################
getNextDecisionEpoch = function(current_decision_epoch) {
  return(
    ifelse(
      current_decision_epoch < 2,
      yes = current_decision_epoch + 0.25,
      no = current_decision_epoch + 0.5
    )
  )
}

#This definition leads to 10% risk threshold as MDP
generateReward = function(current_state, action) {
  if(current_state==AT){
    return(0)
  }else if(current_state==G7){
    return(ifelse(action==BIOPSY, yes = 10, no = -0.2))
  }else{
    return(ifelse(action==BIOPSY, yes = 1, no = 2.1335))
  }
}

#We have to generate data between obs_generation_start_time(not included) 
#and obs_generation_end_time(included)
generateObservation = function(obs_generation_start_time, 
                               obs_generation_end_time, patient_df, 
                               latest_survival_time, earliest_failure_time) {
  
  visitTimeYears_psa = PSA_CHECK_UP_TIME[PSA_CHECK_UP_TIME > obs_generation_start_time & PSA_CHECK_UP_TIME<=obs_generation_end_time]
  visitTimeYears_dre = DRE_CHECK_UP_TIME[DRE_CHECK_UP_TIME > obs_generation_start_time & DRE_CHECK_UP_TIME<=obs_generation_end_time]
  #In the first 2 years this may happen. e.g. generate DRE between 1 and 1.25
  if(length(visitTimeYears_dre)==0){
    return(patient_df[0,])
  }
  
  predicted_outcomes = predictLongitudinalOutcome(fitted_JM, patient_df, 
                                                  latest_survival_time, earliest_failure_time, visitTimeYears_psa, M=0)
  
  n_visitTimeYears_dre = length(visitTimeYears_dre)
  patient_df_return = patient_df[rep(1, n_visitTimeYears_dre),]
  patient_df_return$visitTimeYears = visitTimeYears_dre
  dre_mean_probs = if(n_visitTimeYears_dre==1){
    mean(plogis(predicted_outcomes$predicted_dre_log_odds))
  }else{
    apply(plogis(predicted_outcomes$predicted_dre_log_odds), 1, mean)
  }
  patient_df_return$high_dre = sapply(dre_mean_probs, FUN = rbinom, n=1, size=1)
  
  return(patient_df_return)
}