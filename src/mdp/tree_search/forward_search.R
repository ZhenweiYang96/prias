selectAction = function(patient_df, current_decision_epoch,
                        latest_survival_time, earliest_failure_time,
                        nr_discretized_obs,
                        max_decision_epoch) {
  
  available_actions = if((current_decision_epoch - latest_survival_time)>=MIN_BIOPSY_GAP){
    ACTION_VECTOR
  }else{
    WAIT
  }
  
  #For action wait, the latest survival and earliest failure time 
  #do not change. Thus current belief and future observations can 
  #be sampled from the same distribution. following is a programming trick to
  #speed up our operations.
  
  #First find next decision epoch
  next_decision_epoch = getNextDecisionEpoch(current_decision_epoch)
  
  #The following not only gives us the current belief
  #but also gives us the future observations. addRandomError will get us 
  #samples from the posterior predictive distribution of Y
  expected_future_ifwait = getExpectedFutureOutcomes(fitted_JM, patient_df, 
                                                     latest_survival_time, earliest_failure_time,
                                                     survival_predict_times = current_decision_epoch,
                                                     psa_predict_times = next_decision_epoch,
                                                     dre_predict_times = next_decision_epoch,
                                                     addRandomError = T,
                                                     M=N_MCMC_ITER)
  
  G6_prob = mean(expected_future_ifwait$predicted_surv_prob[1,])
  
  if(BIOPSY %in% available_actions & next_decision_epoch <= max_decision_epoch){
    #The following is the expected future of biopsy when the current state is not G7
    expected_future_ifbiopsy = getExpectedFutureOutcomes(fitted_JM, patient_df, 
                                                         latest_survival_time=current_decision_epoch, 
                                                         earliest_failure_time,
                                                         psa_predict_times = next_decision_epoch,
                                                         dre_predict_times = next_decision_epoch,
                                                         addRandomError = T,
                                                         M=N_MCMC_ITER)
  }
  
  optimal_action = NA
  optimal_reward = -Inf
  
  for(current_action in available_actions){
    current_reward = G6_prob * getReward(G6, current_action) + 
      (1-G6_prob) * getReward(G7, current_action)
    
    if(next_decision_epoch <= max_decision_epoch){
      #probability of not transitioning to AT, given the current action
      non_transition_prob = ifelse(current_action==BIOPSY, G6_prob, 1)
      future_latest_survival_time = ifelse(current_action==BIOPSY, 
                                           current_decision_epoch, 
                                           latest_survival_time)
      
      #Trying to read the MCMC from the back to avoid the burn-in part at the beginning
      for(o in N_MCMC_ITER:(N_MCMC_ITER-nr_discretized_obs + 1)){
        future_data=patient_df[1,]
        future_data$visitTimeYears = next_decision_epoch
        
        future_data$log2psaplus1 = if(current_action==BIOPSY){
          expected_future_ifbiopsy$predicted_psa[1,o]
        }else{
          expected_future_ifwait$predicted_psa[1,o]
        }
        
        future_data$high_dre = if(current_action==BIOPSY){
          rbinom(1,1,expected_future_ifbiopsy$predicted_dre_prob[1,o])
        }else{
          rbinom(1,1,expected_future_ifwait$predicted_dre_prob[1,o])
        }
        
        #This is a conditional one. if action is biopsy, the following is conditional upon not detecting anything currently
        #becuase in the other case, you reach a state AT and future reward is 0
        future_action_reward = selectAction(patient_df=rbind(patient_df, future_data),
                                            current_decision_epoch=next_decision_epoch, 
                                            future_latest_survival_time, earliest_failure_time,
                                            nr_discretized_obs, max_decision_epoch)
        
        current_reward = current_reward + (1/nr_discretized_obs) * non_transition_prob *
          DISCOUNT_FACTORS[as.character(next_decision_epoch)] * future_action_reward$optimal_reward 
        #with a probability 1-non_transition_prob, we get reward 0 because of being in state AT,
        #since that doesn't contribute anything, we do not consider them.
      }
    }
    
    if(current_reward > optimal_reward){
      optimal_action = current_action
      optimal_reward = current_reward
    }
  }
  
  return(list(optimal_action=optimal_action, optimal_reward=optimal_reward))
}