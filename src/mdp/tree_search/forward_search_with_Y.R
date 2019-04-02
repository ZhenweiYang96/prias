selectAction = function(patient_df, current_decision_epoch, 
                        latest_survival_time, earliest_failure_time,
                        max_decision_epoch, cur_biopsies, max_biopsies) {
  
  available_actions = if((current_decision_epoch - latest_survival_time) >= MIN_BIOPSY_GAP){
    ACTION_VECTOR
  }else{
    WAIT
  }

  #First find next decision epoch
  next_decision_epoch = getNextDecisionEpoch(current_decision_epoch)
  expected_future_ifwait = getExpectedFuture_Cat(fitted_JM, patient_df,
                                                   lower_upper_psa_limits = LOWER_UPPER_PSA_LIMITS,
                                                   latest_survival_time, earliest_failure_time,
                                                   survival_predict_times = current_decision_epoch,
                                                   Y_predict_times = next_decision_epoch,
                                                   M=N_MCMC_ITER)
  
  G6_prob = expected_future_ifwait$predicted_surv_prob
  
  if(BIOPSY %in% available_actions & next_decision_epoch <= max_decision_epoch){
    #The following is the expected future of biopsy when the current state is not G7
    expected_future_ifbiopsy = getExpectedFuture_Cat(fitted_JM, patient_df, 
                                                     lower_upper_psa_limits = LOWER_UPPER_PSA_LIMITS,
                                                     latest_survival_time=current_decision_epoch, 
                                                     earliest_failure_time,
                                                     Y_predict_times = next_decision_epoch,
                                                     M=N_MCMC_ITER)
  }
  
  optimal_action = NA
  optimal_reward = -Inf
  
  for(current_action in available_actions){
    current_reward = G6_prob * getReward(G6, current_action, 
                                         current_decision_epoch, latest_survival_time, 
                                         cur_biopsies) +
      (1-G6_prob) * getReward(G7, current_action,
                              current_decision_epoch, latest_survival_time, 
                              cur_biopsies)
    
    fut_optimal_action_chain = vector("list", length(OUTCOME_CAT_NAMES))
    names(fut_optimal_action_chain) = OUTCOME_CAT_NAMES
    
    if(next_decision_epoch <= max_decision_epoch){
      #probability of not transitioning to AT, given the current action
      non_AT_transition_prob = ifelse(current_action==BIOPSY, G6_prob, 1)
      future_latest_survival_time = ifelse(current_action==BIOPSY, 
                                           current_decision_epoch, 
                                           latest_survival_time)
      
      for(o in 1:length(OUTCOME_CAT_NAMES)){
        Y_prob = ifelse(current_action==BIOPSY,
                        expected_future_ifbiopsy$predicted_Y_prob[o,1],
                        expected_future_ifwait$predicted_Y_prob[o,1])
        
        if(Y_prob > 0){
          future_data = patient_df[1,]
          future_data$visitTimeYears = next_decision_epoch
          future_data$psa_cat_data = T
          future_data[, c("log2psaplus1", "high_dre")] = OUTCOME_PSA_DRE_CAT[o, ]
          
          future_cur_biopsies = ifelse(current_action==BIOPSY, cur_biopsies + 1, cur_biopsies)
          future_action_reward = selectAction(patient_df=rbind(patient_df, future_data),
                                              current_decision_epoch=next_decision_epoch,
                                              future_latest_survival_time, earliest_failure_time,
                                              max_decision_epoch, future_cur_biopsies, max_biopsies)
          
          fut_optimal_action_chain[[OUTCOME_CAT_NAMES[o]]] = future_action_reward$optimal_action_chain
          
          current_reward = current_reward + DISCOUNT_FACTORS[as.character(next_decision_epoch)] *
              non_AT_transition_prob * Y_prob * future_action_reward$optimal_reward
          #with a probability 1-non_AT_transition_prob, we get reward 0 because of being in state AT,
          #since that doesn't contribute anything, we do not consider them.
        }
      }
    }
    
    if(current_reward > optimal_reward){
      optimal_action = current_action
      
      if(all(sapply(fut_optimal_action_chain, is.null))){
        optimal_action_chain = list(current_action)
        names(optimal_action_chain) = "Now"
      }else{
        optimal_action_chain = list(current_action, fut_optimal_action_chain)
        names(optimal_action_chain) = c("Now", "Future")
      }
      optimal_reward = current_reward
    }
  }
  
  return(list(optimal_action=optimal_action, 
              optimal_reward=optimal_reward,
              optimal_action_chain = optimal_action_chain))
}
