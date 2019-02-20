selectAction = function(patient_df, current_decision_epoch, G6_probs=NULL,
                        latest_survival_time, earliest_failure_time,
                        max_decision_epoch) {
  
  available_actions = if((current_decision_epoch - latest_survival_time)>=MIN_BIOPSY_GAP){
    ACTION_VECTOR
  }else{
    WAIT
  }
  
  if(is.null(G6_probs)){
    survival_predict_times = PSA_CHECK_UP_TIME[PSA_CHECK_UP_TIME>=current_decision_epoch & PSA_CHECK_UP_TIME<=max_decision_epoch]
    #The following not only gives us the current belief
    cur_state = getExpectedFutureOutcomes(fitted_JM, patient_df, 
                                          latest_survival_time, earliest_failure_time,
                                          survival_predict_times,
                                          M=N_MCMC_ITER)
    
    G6_probs = apply(cur_state$predicted_surv_prob, 1, mean)
  }
  
  G6_prob = G6_probs[1]
  G6_probs = G6_probs[-1]
  
  optimal_action = NA
  optimal_reward = -Inf
  
  next_decision_epoch = getNextDecisionEpoch(current_decision_epoch)
  for(current_action in available_actions){
    current_reward = G6_prob * getReward(G6, current_action) +
        (1-G6_prob) * getReward(G7, current_action)
    # if(current_action==BIOPSY){
    #    current_reward = intercept_b + (1-G6_prob) * slope_b
    # }else{
    #    current_reward = intercept_w + (1-G6_prob) * slope_w
    # }
    
    if(next_decision_epoch <= max_decision_epoch){
      #probability of not transitioning to AT, given the current action
      non_transition_prob = ifelse(current_action==BIOPSY, G6_prob, 1)
      future_latest_survival_time = ifelse(current_action==BIOPSY, 
                                           current_decision_epoch, 
                                           latest_survival_time)
      
      future_G6_probs = if(current_action==BIOPSY){
        NULL
      }else{
        G6_probs
      }
      #This is a conditional one. if action is biopsy, the following is conditional upon not detecting anything currently
      #becuase in the other case, you reach a state AT and future reward is 0
      future_action_reward = selectAction(patient_df, next_decision_epoch, future_G6_probs,
                                          future_latest_survival_time, earliest_failure_time,
                                          max_decision_epoch)
      
      current_reward = current_reward + DISCOUNT_FACTORS[as.character(next_decision_epoch)] * 
        non_transition_prob * future_action_reward$optimal_reward 
      #with a probability 1-non_transition_prob, we get reward 0 because of being in state AT,
      #since that doesn't contribute anything, we do not consider them.
    }
    
    if(current_reward > optimal_reward){
      optimal_action = current_action
      optimal_reward = current_reward
    }
  }
  
  return(list(optimal_action=optimal_action, optimal_reward=optimal_reward,
              G6_probs=c(G6_prob,G6_probs)))
}
