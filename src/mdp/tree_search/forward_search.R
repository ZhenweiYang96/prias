generateFutureObs = function(patient_df, belief, action, current_decision_epoch, 
                             latest_survival_time, earliest_failure_time) {
  future_decision_epoch = getNextDecisionEpoch(current_decision_epoch, action)
  
  if(action==BIOPSY){
    #If current action is BIOPSY then 
    #probability of NA at the future decision epoch is equal to the current belief in G7
    if(rbinom(1,1,prob = belief[G7])==1){
      future_data = patient_df[0,]
      future_data_prob = belief[G7]
    }else{
      future_data_with_prob = generateObservation(current_decision_epoch, future_decision_epoch,
                                                  patient_df, current_decision_epoch, earliest_failure_time)
      future_data = future_data_with_prob$future_data
      future_data_prob = belief[G6] * future_data_with_prob$future_data_prob
    }
  }else{
    future_data_with_prob = generateObservation(current_decision_epoch, future_decision_epoch,
                                                patient_df, latest_survival_time, earliest_failure_time)
    future_data = future_data_with_prob$future_data
    future_data_prob = future_data_with_prob$future_data_prob
  }
  
  return(list(future_decision_epoch=future_decision_epoch,
              future_data=future_data, future_data_prob=future_data_prob))
}

selectAction = function(patient_df, current_decision_epoch,
                        latest_survival_time, earliest_failure_time,
                        discount_factors,
                        nr_discretized_obs,
                        max_decision_epoch) {
  
  ##Stopping rule
  reached_AT = is.na(tail(patient_df$log2psaplus1, 1))
  if(current_decision_epoch > max_decision_epoch | reached_AT){
    return(list(optimal_action = NA, optimal_reward=0))
  }else{
    G6_prob = survfitJM_noSeed(fitted_JM, 
                               patient_df, idVar = "P_ID", 
                               last.time = latest_survival_time,
                               survTimes = current_decision_epoch)[[1]][[1]][,"Mean"]
    current_belief = c(G6=G6_prob, G7 = 1-G6_prob, AT = 0)
    
    optimal_action = NA
    optimal_reward = -Inf
    
    for(current_action in ACTION_VECTOR){
      current_reward = sum(current_belief * sapply(names(current_belief), FUN = generateReward, action=current_action))
      
      for(o in 1:nr_discretized_obs){
        future_obs = generateFutureObs(patient_df, current_belief, current_action, 
                                       current_decision_epoch, latest_survival_time, earliest_failure_time)
        
        future_action_reward = selectAction(rbind(patient_df, future_obs$future_data),
                                                 future_obs$future_decision_epoch, 
                                                 latest_survival_time, earliest_failure_time,
                                                 discount_factors,
                                                 nr_discretized_obs, max_decision_epoch)
        current_reward = current_reward + 
          discount_factors[as.character(future_obs$future_decision_epoch)] * future_obs$future_data_prob * future_action_reward$reward
      }
      
      if(current_reward > optimal_reward){
        optimal_action = current_action
        optimal_reward = current_reward
      }
    }
    
    return(list(optimal_action=optimal_action, optimal_reward=optimal_reward))
  }
}