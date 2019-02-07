selectAction = function(patient_df, current_decision_epoch,
                        latest_survival_time, earliest_failure_time,
                        nr_discretized_obs,
                        max_decision_epoch) {
  
  available_actions = ifelse((current_decision_epoch - earliest_failure_time)>=MIN_BIOPSY_GAP, ACTION_VECTOR, WAIT)
  
  #For action wait, the latest survival and earliest failure time 
  #do not change. Thus current belief and future observations can 
  #be sampled from the same distribution. following is a programming trick to
  #speed up our operations.
  
  #First find next decision epoch
  next_decision_epoch = getNextDecisionEpoch(current_decision_epoch)
  
  #The following not only gives us the current belief
  #but also gives us the future observations. ignoreErrorTerm will get us 
  #samples from the posterior predictive distribution of Y
  expected_future_wait = getExpectedFutureOutcomes(fitted_JM, patient_df, 
                                                   latest_survival_time, earliest_failure_time,
                                                   survival_predict_times = current_decision_epoch,
                                                   psa_predict_times = next_decision-epoch,
                                                   dre_predict_times = next_decision-epoch,
                                                   ignoreErrorTerm = F,
                                                   M=N_MCMC_ITER)
  
  G6_prob = mean(expected_future_wait$predicted_surv_prob[1,])
  current_belief = c(G6=G6_prob, G7 = 1-G6_prob, AT = 0)
  
  if(BIOPSY %in% available_actions & next_decision_epoch >= max_decision_epoch){
    #We also do a shortcut for future obs when action is BIOPSY
    #The following is the expected future of biopsy when current state is not G7
    expected_future_biopsy = getExpectedFutureOutcomes(fitted_JM, patient_df, 
                                                       latest_survival_time=current_decision_epoch, 
                                                       earliest_failure_time,
                                                       psa_predict_times = next_decision_epoch,
                                                       dre_predict_times = next_decision_epoch,
                                                       ignoreErrorTerm = F,
                                                       M=N_MCMC_ITER)
  }
  
  optimal_action = NA
  optimal_reward = -Inf
  
  for(current_action in available_actions){
    current_reward = sum(current_belief * sapply(names(current_belief), FUN = generateReward, action=current_action))
    
    if(next_decision_epoch >= max_decision_epoch){
      #probability of not transitioning to AT, given the current action
      non_transition_prob = ifelse(current_action==BIOPSY, G6_prob, 1)
      
      for(o in 1:nr_discretized_obs){
        future = list(future_decision_epoch=next_decision_epoch,
                      future_data=patient_df[1,])
        future$future_data$visitTimeYears = next_decision_epoch
        
        #Trying to read the MCMC from the back to avoid the burn-in part at the beginning
        obs_index = N_MCMC_ITER - o - 1
        if(current_action==BIOPSY){
          future$future_data$log2psaplus1 = expected_future_biopsy$predicted_psa[1,obs_index]
          future$future_data$high_dre = rbinom(1,1,expected_future_biopsy$predicted_dre_prob[1,obs_index])
          
        }else{
          future$future_data$log2psaplus1 = expected_future_wait$predicted_psa[1,obs_index]
          future$future_data$high_dre = rbinom(1,1,expected_future_wait$predicted_dre_prob[1,obs_index])
        }
        
        future_latest_survival_time = ifelse(current_action==BIOPSY, current_decision_epoch, latest_survival_time)
        future_action_reward = selectAction(patient_df=rbind(patient_df, future$future_data),
                                            current_decision_epoch=next_decision_epoch, 
                                            future_latest_survival_time, earliest_failure_time,
                                            nr_discretized_obs, max_decision_epoch)
        
        current_reward = current_reward + (1/nr_discretized_obs) * non_transition_prob *
          DISCOUNT_FACTORS[as.character(future$future_decision_epoch)] * future_action_reward$optimal_reward 
        
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