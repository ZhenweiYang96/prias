selectAction = function(patient_df, current_decision_epoch, G6_probs=NULL,
                        latest_survival_time, earliest_failure_time,
                        max_decision_epoch, cur_biopsies=0, max_biopsies=Inf) {
  available_actions = if((current_decision_epoch - latest_survival_time)>=MIN_BIOPSY_GAP & cur_biopsies<max_biopsies){
    ACTION_VECTOR
  }else{
    WAIT
  }
  
  if(is.null(G6_probs)){
    survival_predict_times = BIOPSY_TEST_TIMES[BIOPSY_TEST_TIMES>=current_decision_epoch & BIOPSY_TEST_TIMES<=max_decision_epoch]
    
    cur_state = getExpectedFutureOutcomes(fitted_JM, patient_df, 
                                          latest_survival_time, earliest_failure_time,
                                          survival_predict_times,
                                          M=N_MCMC_ITER)
    
    G6_probs = apply(cur_state$predicted_surv_prob, 1, mean)
  }
  
  G6_prob = G6_probs[1]
  
  optimal_action = NA
  optimal_reward = -Inf
  
  next_decision_epoch = getNextDecisionEpoch(current_decision_epoch)
  for(current_action in available_actions){
    current_reward = G6_prob * getReward(G6, current_action, 
                                         current_decision_epoch, latest_survival_time, 
                                         cur_biopsies) +
      (1-G6_prob) * getReward(G7, current_action, 
                              current_decision_epoch, latest_survival_time, 
                              cur_biopsies)
    
    fut_optimal_action_chain = ""
    if(next_decision_epoch <= max_decision_epoch){
      #probability of not transitioning to AT, given the current action
      non_transition_prob = ifelse(current_action==BIOPSY, G6_prob, 1)
      future_latest_survival_time = ifelse(current_action==BIOPSY, 
                                           current_decision_epoch, 
                                           latest_survival_time)
      
      future_G6_probs = if(current_action==BIOPSY){
        NULL
      }else{
        G6_probs[-1]
      }
      
      future_cur_biopsies = ifelse(current_action==BIOPSY, cur_biopsies + 1, cur_biopsies)
      #This is a conditional one. if action is biopsy, the following is conditional upon not detecting anything currently
      #becuase in the other case, you reach a state AT and future reward is 0
      future_action_reward = selectAction(patient_df, next_decision_epoch, future_G6_probs,
                                          future_latest_survival_time, earliest_failure_time,
                                          max_decision_epoch, future_cur_biopsies, max_biopsies)
      
      fut_optimal_action_chain = future_action_reward$optimal_action_chain
      
      current_reward = current_reward + DISCOUNT_FACTORS[as.character(next_decision_epoch)] * 
        non_transition_prob * future_action_reward$optimal_reward 
      #with a probability 1-non_transition_prob, we get reward 0 because of being in state AT,
      #since that doesn't contribute anything, we do not consider them.
    }
    
    if(current_reward > optimal_reward){
      optimal_action = current_action
      optimal_action_chain = paste(current_action, fut_optimal_action_chain)
      optimal_reward = current_reward
    }
  }
  
  # print(paste("Latest survival time:",latest_survival_time,
  #             "Current Time:", current_decision_epoch, 
  #             "--- Optimal action:", optimal_action,
  #             "--- Optimal reward:", optimal_reward))
  
  return(list(optimal_action=optimal_action, optimal_action_chain=optimal_action_chain,
              optimal_reward=optimal_reward,
              G6_probs=G6_probs))
}
