#Constants
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

EXPLORED_TREE = list()

#CONSTANT STRINGS
VISIT_COUNT_BIOPSY = "VCB"
Q_VALUE_BIOPSY = "QVB"
VISIT_COUNT_WAIT = "VCW"
Q_VALUE_WAIT = "QVW"
G6_COUNT = "G6C"

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
        yes = current_decision_epoch + 0.5,
        no = current_decision_epoch + 0.5
      )
    )
  }
}

#This definition leads to 10% risk threshold as MDP
generateReward = function(action, current_state, current_decision_epoch,
                          latest_survival_time, earliest_failure_time) {
  if(current_state==AT){
    return(0)
  }else if(current_state==G7){
    return(ifelse(action==BIOPSY, yes = 10, no = -0.2))
  }else{
    return(ifelse(action==BIOPSY, yes = 1, no = 2.1335))
  }
}

# generateReward = function(action, current_state, current_decision_epoch, 
#                           latest_survival_time, earliest_failure_time) {
#   if(current_state==AT){
#     return(0)
#   }else if(current_state==G7){
#     estimated_unknown_delay = 0.5 * (earliest_failure_time - latest_survival_time)
#     next_decision_epoch_wait = getNextDecisionEpoch(current_decision_epoch, WAIT)
#     
#     #THe no part could be weighted, using Infinite series formula??
#     return(ifelse(action==BIOPSY, 
#                   yes = 1/(estimated_unknown_delay + current_decision_epoch - earliest_failure_time), 
#                   no = 1/(estimated_unknown_delay + current_decision_epoch + next_decision_epoch_wait - earliest_failure_time))) 
#   }else{
#     return(ifelse(action==BIOPSY, yes = 0.1, no = 0.21335))
#   }
# }

generateObservation = function(current_state, obs_generation_start_time, 
                               obs_generation_end_time, patient_obs_df, 
                               latest_survival_time, earliest_failure_time) {
  
  # return(patient_obs_df[0,])
  
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
  #In the first 2 years this may happen. e.g. generate DRE between 1 and 1.25
  if(length(visitTimeYears_dre)==0){
    return(patient_obs_df[0,])
  }
  
  #currently only for DRE
  predicted_outcomes = predictLongitudinalOutcome(fitted_JM, patient_obs_df, 
                                                  latest_survival_time, earliest_failure_time, visitTimeYears_dre, M=0)
  
  n_visitTimeYears_dre = length(visitTimeYears_dre)
  patient_obs_df_return = patient_obs_df[rep(1, n_visitTimeYears_dre),]
  patient_obs_df_return$visitTimeYears = visitTimeYears_dre
  dre_mean_probs = if(n_visitTimeYears_dre==1){
    mean(plogis(predicted_outcomes$predicted_dre_log_odds))
  }else{
    apply(plogis(predicted_outcomes$predicted_dre_log_odds), 1, mean)
  }
  patient_obs_df_return$high_dre = sapply(dre_mean_probs, FUN = rbinom, n=1, size=1)
  
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
    futureEpochSurvivalProb = survfitJM_noSeed(fitted_JM, 
                                               patient_obs_df, idVar = "P_ID", 
                                               last.time = latest_survival_time,
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
                   latest_survival_time, earliest_failure_time, max_decision_epoch, discount_factors) {
  
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
  
  action_immediate_reward = generateReward(action, current_state, 
                                           current_decision_epoch, 
                                           latest_survival_time, earliest_failure_time)
  
  rollout_reward = rollout(future$future_state,
                           future$future_decision_epoch,
                           future$future_patient_obs_df,
                           latest_survival_time,
                           earliest_failure_time,
                           max_decision_epoch, discount_factors)
  
  return(action_immediate_reward + discount_factors[as.character(current_decision_epoch)] * rollout_reward)
}


simulateMCTS = function(current_state,
                        current_decision_epoch,
                        patient_obs_df,
                        latest_survival_time,
                        earliest_failure_time,
                        node_hash,
                        max_decision_epoch, discount_factors) {
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
    EXPLORED_TREE[[node_hash]] <<- c(rep(0,4), ifelse(current_state==G6,1,0))
    names(EXPLORED_TREE[[node_hash]]) <<- c(VISIT_COUNT_BIOPSY, Q_VALUE_BIOPSY, 
                                            VISIT_COUNT_WAIT, Q_VALUE_WAIT, G6_COUNT)
    return(rollout(current_state, current_decision_epoch, patient_obs_df, 
                   latest_survival_time, earliest_failure_time, max_decision_epoch, discount_factors))
  }
  
  current_visit_count_biopsy = EXPLORED_TREE[[node_hash]][VISIT_COUNT_BIOPSY]
  current_q_value_biopsy = EXPLORED_TREE[[node_hash]][Q_VALUE_BIOPSY]
  current_visit_count_wait = EXPLORED_TREE[[node_hash]][VISIT_COUNT_WAIT]
  current_q_value_wait = EXPLORED_TREE[[node_hash]][Q_VALUE_WAIT]
  
  EXPLORED_TREE[[node_hash]][G6_COUNT] <<- EXPLORED_TREE[[node_hash]][G6_COUNT] + ifelse(current_state==G6,1,0)
  
  #GO GREEDY 1
  #This algorithm can be reduced to the one below GO GREEDY 2
  #total_current_visit_count = current_visit_count_biopsy + current_visit_count_wait
  #reward_action_biopsy = current_q_value_biopsy + EXPLORATION_REWARD_CONSTANT * sqrt(log(total_current_visit_count)/current_visit_count_biopsy)
  #reward_action_wait = current_q_value_wait + EXPLORATION_REWARD_CONSTANT * sqrt(log(total_current_visit_count)/current_visit_count_wait)
  
  #GO GREEDY 2
  #reward_action_biopsy = current_q_value_biopsy + EXPLORATION_REWARD_CONSTANT * sqrt(1/current_visit_count_biopsy)
  #reward_action_wait = current_q_value_wait + EXPLORATION_REWARD_CONSTANT * sqrt(1/current_visit_count_wait)
  #action = ifelse(reward_action_biopsy <= reward_action_wait, yes = WAIT, no = BIOPSY)
  
  #RANDOM ACTION
  action = sample(c(BIOPSY, WAIT), size = 1)
  future = generateFuture(current_state, current_decision_epoch,
                          action, patient_obs_df, 
                          latest_survival_time, earliest_failure_time)
  
  action_immediate_reward = generateReward(action, current_state,
                                           current_decision_epoch, 
                                           latest_survival_time, earliest_failure_time)
  
  new_node_hash = paste(node_hash, future$future_node_hash)
  #new_node_hash = paste(node_hash, action)
  action_q_value = action_immediate_reward + discount_factors[as.character(current_decision_epoch)] * simulateMCTS(future$future_state,
                                                                            future$future_decision_epoch,
                                                                            future$future_patient_obs_df,
                                                                            latest_survival_time, earliest_failure_time,
                                                                            new_node_hash,
                                                                            max_decision_epoch, discount_factors)
  if(action==BIOPSY){
    EXPLORED_TREE[[node_hash]][VISIT_COUNT_BIOPSY] <<- current_visit_count_biopsy + 1
    EXPLORED_TREE[[node_hash]][Q_VALUE_BIOPSY] <<- (action_q_value + current_visit_count_biopsy * current_q_value_biopsy)/(current_visit_count_biopsy + 1)
  }else{
    EXPLORED_TREE[[node_hash]][VISIT_COUNT_WAIT] <<- current_visit_count_wait + 1
    EXPLORED_TREE[[node_hash]][Q_VALUE_WAIT] <<- (action_q_value + current_visit_count_wait * current_q_value_wait)/(current_visit_count_wait + 1)
  }
  
  return(action_q_value)
}

selectAction = function(patient_obs_df, n_simulations = 100,
                        current_decision_epoch = 1,
                        latest_survival_time = 0,
                        max_decision_epoch = MAX_FOLLOW_UP_TIME) {
  
  decision_epochs = PSA_CHECK_UP_TIME[PSA_CHECK_UP_TIME>=current_decision_epoch & PSA_CHECK_UP_TIME<=max_decision_epoch]
  discount_factors = DISCOUNT_FACTOR^(1:length(decision_epochs))
  names(discount_factors) = decision_epochs
  
  EMPTY_HISTORY = "-"
  #TREE expressed as a LIST
  EXPLORED_TREE <<- list()
  G6_prob = survfitJM_noSeed(fitted_JM, 
                             patient_obs_df, idVar = "P_ID", 
                             last.time = latest_survival_time,
                             survTimes = current_decision_epoch)[[1]][[1]][,"Mean"]
  
  print(discount_factors)
  G6_prob = 0.7
  starting_states = sample(c(G6, G7), size = n_simulations, prob = c(G6_prob, 1-G6_prob), replace = T)
  
  for (i in 1:n_simulations) {
    simulateMCTS(
      starting_states[i], current_decision_epoch, patient_obs_df,
      latest_survival_time, earliest_failure_time=Inf,
      EMPTY_HISTORY,
      max_decision_epoch, discount_factors
    )
  }
  
  final_action = ifelse(EXPLORED_TREE[[EMPTY_HISTORY]][Q_VALUE_BIOPSY] > EXPLORED_TREE[[EMPTY_HISTORY]][Q_VALUE_WAIT], BIOPSY, WAIT)
  explored_tree = EXPLORED_TREE
  EXPLORED_TREE <<- list()
  return(list(final_action=final_action, explored_tree=explored_tree))
}