fitted_JM = NULL
patient_data_starting = NULL

#Maximum follow up time in years
MAX_FOLLOW_UP_TIME = 10

#Minimum biopsy gap should not be less than 0.25 years
MIN_BIOPSY_GAP = 1

#Decision Epochs in years
PSA_CHECK_UP_TIME = c(seq(0, 2, 0.25), seq(2.5, MAX_FOLLOW_UP_TIME, 0.5))

#DRE check up time years
DRE_CHECK_UP_TIME = c(seq(0, MAX_FOLLOW_UP_TIME, 0.5))

#Actions
BIOPSY = "B"
WAIT = "W"
ACTION_VECTOR = c(BIOPSY, WAIT)

#States
G7 = "G7"
G6 = "G6"
AT = "AT"

EMPTY_HISTORY = "-"

#TREE expressed as a LIST
EXPLORED_TREE = list()
#CONSTANT STRINGS
VISIT_COUNT_BIOPSY = "VCB"
Q_VALUE_BIOPSY = "QVB"
VISIT_COUNT_WAIT = "VCW"
Q_VALUE_WAIT = "QVW"

#Must be chosen suitably
EXPLORATION_REWARD_CONSTANT = 10
DISCOUNT_FACTOR = 1

#########################################
# We define various functions from here onwards
#########################################

getNextDecisionEpoch = function(current_decision_epoch, current_action) {
  if (current_action == BIOPSY) {
    return(current_decision_epoch + MIN_BIOPSY_GAP)
  } else{
    return(
      ifelse(
        current_decision_epoch < 2,
        yes = current_decision_epoch + 0.25,
        no = current_decision_epoch + 0.5
      )
    )
  }
}

generateReward = function(action, action_state) {
  if(action_state==AT){
    return(0)
  }else if(action_state==G7){
    #10 points for correct biopsy, -10 for unnecessary waiting
    return(ifelse(action==BIOPSY, yes = 10, no = -10))
  }else{
    #10 points for correct waiting, -10 for unnecessary biopsy
    return(ifelse(action==BIOPSY, yes = -10, no = 10))
  }
}

generateObservation = function(current_state, obs_generation_start_time, 
                               obs_generation_end_time, patient_obs_df, 
                               latest_survival_time, earliest_failure_time) {
  if (current_state == AT) {
    #Return empty row that has no effect on rbind
    return(patient_obs_df[0,])
  } 
  
  #OTHERWISE GENERATE DATA AT THE CURRENT EPOCH AS PER
  #Obs patient data, latest_survival_time, and earliest_failure_time
  #We have to generate data between obs_generation_start_time(not included) 
  #and obs_generation_end_time(included)
  
  #currently only for DRE
  #visitTimeYears_psa = PSA_CHECK_UP_TIME[PSA_CHECK_UP_TIME > obs_generation_start_time & PSA_CHECK_UP_TIME<=obs_generation_end_time]
  visitTimeYears_dre = DRE_CHECK_UP_TIME[DRE_CHECK_UP_TIME > obs_generation_start_time & DRE_CHECK_UP_TIME<=obs_generation_end_time]
  
  #currently only for DRE
  predicted_outcomes = predictLongitudinalOutcome(fitted_JM, patient_obs_df, 
                                                  latest_survival_time, earliest_failure_time, visitTimeYears_dre)
  
  patient_obs_df_return = patient_obs_df[rep(1, length(visitTimeYears_dre)),]
  patient_obs_df_return$visitTimeYears = visitTimeYears_dre
  patient_obs_df_return$high_dre = sapply(plogis(predicted_outcomes$predicted_dre_log_odds), FUN = sample, x=c(1,0), size=1, replace=T)
  
  return(patient_obs_df_return)
}

#Return future state, a new observation, and reward of action given current state
#latest_survival_time means the last time the state was G6. 
#This is used to generate patient data when current state is G7 and we need to know the interval censoring time
generateFuture = function(current_state, current_decision_epoch, action,
                          patient_obs_df, latest_survival_time, 
                          earliest_failure_time) {
  future_decision_epoch = getNextDecisionEpoch(current_decision_epoch, action)
  
  if (current_state == AT) {
    future_state = AT
  } else if (current_state == G7){
    future_state = ifelse(action == BIOPSY, yes = AT, no = G7)
  } else {
    #The current state is G6 
    #The current_decision_epoch and latest_survival_time are equal
    futureEpochSurvivalProb = survfitJM(fitted_JM, 
                                        patient_obs_df, idVar = "P_ID", 
                                        last.time = current_decision_epoch,
                                        survTimes = future_decision_epoch)[[1]][[1]][,"Mean"]
    
    future_state = sample(
      c(G6, G7),
      size = 1,
      prob = c(futureEpochSurvivalProb, 1 - futureEpochSurvivalProb)
    )
  }
  
  future_data = generateObservation(future_state, current_decision_epoch,
                                    future_decision_epoch, patient_obs_df,
                                    latest_survival_time, earliest_failure_time)
  
  #add hash of PSA later
  future_node_hash = paste0(action, paste0(future_data$high_dre, collapse = ''))
  
  return(list(
    future_state = future_state,
    future_decision_epoch = future_decision_epoch,
    future_patient_obs_df = rbind(patient_obs_df, future_data),
    future_node_hash = future_node_hash
  ))
}

#Randomly exploring the tree
rollout = function(current_state, current_decision_epoch, patient_obs_df,
                   latest_survival_time, earliest_failure_time, max_decision_epoch) {
  
  #STOPPING RULE
  if (current_decision_epoch > max_decision_epoch | current_state==AT) {
    return(0)
  } 
  
  #OTHERWISE THE FOLLOWING EXECUTES
  if(current_state==G6) {
    latest_survival_time = current_decision_epoch
  }else if(is.infinite(earliest_failure_time)){
    #This condition will be trigerred only once, when patient enters G7 state for the first time
    earliest_failure_time = current_decision_epoch
  }
  
  #Take random action or according to function you desire. I prefer random
  action = sample(c(BIOPSY, WAIT), size = 1)
  
  future = generateFuture(current_state, current_decision_epoch,
                          action, patient_obs_df, 
                          latest_survival_time, earliest_failure_time)
  
  action_immediate_reward = generateReward(action, current_state)
  
  rollout_reward = rollout(future$future_state,
                           future$future_decision_epoch,
                           future$future_patient_obs_df,
                           latest_survival_time,
                           earliest_failure_time,
                           max_decision_epoch)
  
  return(action_immediate_reward + DISCOUNT_FACTOR * rollout_reward)
}


simulateMCTS = function(current_state,
                        current_decision_epoch,
                        patient_obs_df,
                        latest_survival_time,
                        earliest_failure_time,
                        node_hash,
                        max_decision_epoch) {
  #STOPPING RULE
  if (current_decision_epoch > max_decision_epoch | current_state==AT) {
    return(0)
  } 
  
  #OTHERWISE THE FOLLOWING EXECUTES
  if(current_state==G6) {
    latest_survival_time = current_decision_epoch
  }else if(is.infinite(earliest_failure_time)){
    #This condition will be trigerred only once, when patient enters G7 state for the first time
    earliest_failure_time = current_decision_epoch
  }
  
  if(is.null(EXPLORED_TREE[[node_hash]])){
    #GO RANDOM USING ROLLOUTS
    EXPLORED_TREE[[node_hash]] = c(VISIT_COUNT_BIOPSY=0, Q_VALUE_BIOPSY=0,
                                   VISIT_COUNT_WAIT=0, Q_VALUE_WAIT=0)
    return(rollout(current_state, current_decision_epoch, patient_obs_df, 
                   latest_survival_time, earliest_failure_time, max_decision_epoch))
  }
  
  current_visit_count_biopsy = EXPLORED_TREE[[node_hash]][VISIT_COUNT_BIOPSY]
  current_q_value_biopsy = EXPLORED_TREE[[node_hash]][Q_VALUE_BIOPSY]
  current_visit_count_wait = EXPLORED_TREE[[node_hash]][VISIT_COUNT_WAIT]
  current_q_value_wait = EXPLORED_TREE[[node_hash]][Q_VALUE_WAIT]
  total_current_visit_count = current_visit_count_biopsy + current_visit_count_wait
  
  #GO GREEDY
  reward_action_biopsy = current_q_value_biopsy + EXPLORATION_REWARD_CONSTANT * sqrt(log(total_current_visit_count)/current_visit_count_biopsy)
  reward_action_wait = current_q_value_wait + EXPLORATION_REWARD_CONSTANT * sqrt(log(total_current_visit_count)/current_visit_count_wait)
  
  action = ifelse(reward_action_biopsy > reward_action_wait, yes = BIOPSY, no = WAIT)
  
  future = generateFuture(current_state, current_decision_epoch,
                          action, patient_obs_df, 
                          latest_survival_time, earliest_failure_time)
  
  action_immediate_reward = generateReward(action, current_state)
  
  new_node_hash = paste(node_hash, future$future_node_hash)
  action_q_value = action_immediate_reward + DISCOUNT_FACTOR * simulateMCTS(future$future_state,
                                                                            future$future_decision_epoch,
                                                                            future$future_patient_obs_df,
                                                                            latest_survival_time, earliest_failure_time,
                                                                            new_node_hash,
                                                                            max_decision_epoch)
  
  if(action==BIOPSY){
    EXPLORED_TREE[[node_hash]][VISIT_COUNT_BIOPSY] = current_visit_count_biopsy + 1
    EXPLORED_TREE[[node_hash]][Q_VALUE_BIOPSY] = (action_q_value + current_visit_count_biopsy * current_q_value_biopsy)/(current_visit_count_biopsy + 1)
  }else{
    EXPLORED_TREE[[node_hash]][VISIT_COUNT_WAIT] = current_visit_count_wait + 1
    EXPLORED_TREE[[node_hash]][Q_VALUE_BIOPSY] = (action_q_value + current_visit_count_wait * current_q_value_wait)/(current_visit_count_wait + 1)
  }
  
  return(action_q_value)
}

selectAction = function(patient_obs_df, n_simulations = 100,
                        current_decision_epoch = 1,
                        max_decision_epoch = MAX_FOLLOW_UP_TIME) {
  latest_survival_time = max(patient_obs_df$visitTimeYears[!is.na(patient_obs_df$gleason)])
  
  G6_prob = survfitJM(fitted_JM, 
                      patient_obs_df, idVar = "P_ID", 
                      last.time = latest_survival_time,
                      survTimes = current_decision_epoch)[[1]][[1]][,"Mean"]
  
  for (i in 1:n_simulations) {
    starting_state = sample(c(G6, G7), 1, prob = c(G6_prob, 1-G6_prob))
    simulateMCTS(
      starting_state, current_decision_epoch, patient_obs_df,
      latest_survival_time, earliest_failure_time=Inf,
      EMPTY_HISTORY,
      max_decision_epoch
    )
  }
  
  if (EXPLORED_TREE[[EMPTY_HISTORY]][Q_VALUE_BIOPSY] > EXPLORED_TREE[[EMPTY_HISTORY]][Q_VALUE_WAIT]) {
    return(BIOPSY)
  } else{
    return(WAIT)
  }
}