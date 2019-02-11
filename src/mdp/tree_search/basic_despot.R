# generateFuture = function(state, action, patient_df, current_decision_epoch,
#                           latest_survival_time, earliest_failure_time){
#   next_decision_epoch = getNextDecisionEpoch(current_decision_epoch, action)
#   
#   if(current_state == G7 & action==BIOPSY){
#     future_obs = patient_df[0,]
#   }else{
#     future_obs = generateObservation(current_decision_epoch, next_decision_epoch,
#                                      patient_df, latest_survival_time, 
#                                      earliest_failure_time)
#   }
#   patient_df = rbind(patient_df, future_obs)
#   
#   if(current_state == G7){
#     future_state = ifelse(action==BIOPSY, AT, G7)
#   }else{
#     future_belief = getExpectedFutureOutcomes(fitted_JM, patient_df,
#                                               latest_survival_time, earliest_failure_time,
#                                               survival_predict_times = next_decision_epoch,
#                                               M=N_MCMC_ITER)
#     
#     G6_prob = mean(future_belief$predicted_surv_prob[1,])
#     
#     future_state = sample(x=STATE_VECTOR, size=1, prob=c(G6_prob, 1-G6_prob, 0))
#   }
#   
#   return(list(state = future_state, patient_df = patient_df,
#               decision_epoch = next_decision_epoch))
# }

createNode = function(patient_df, current_decision_epoch, 
                      state_distribution, latest_biopsy_time,
                      latest_survival_time, earliest_failure_time,
                      nr_discretized_obs, max_decision_epoch){
  
  biopsy_children = wait_children = NULL
  
  #Do not check for current_decision_epoch <= max_decision_epoch here
  if(current_decision_epoch < max_decision_epoch){
    next_decision_epoch = getNextDecisionEpoch(current_decision_epoch)
    biopsy_children = list()
    wait_children = list()
    future_outcome_cur_G6 = NULL
    
    if(current_decision_epoch - latest_biopsy_time >= 1){
      #For action Biopsy if current state is G6
      if(state_distribution[G6] > 0){
        future_outcome_cur_G6 = getExpectedFutureOutcomes(fitted_JM, patient_df, 
                                                          latest_survival_time = current_decision_epoch,
                                                          earliest_failure_time = earliest_failure_time,
                                                          dre_predict_times = next_decision_epoch,
                                                          psa_predict_times = next_decision_epoch,
                                                          addRandomError = T, M = N_MCMC_ITER)
        
        #At the moment I resample from the last nr_discretized_obs
        #assuming that nr_discretized_obs represent the distribution of Y well
        obs_index_sample_G6_biopsy = sample(N_MCMC_ITER:(N_MCMC_ITER-nr_discretized_obs + 1), 
                                          size=state_distribution[G6], replace = T)
        obs_distribution_G6_biopsy = table(obs_index_sample_G6_biopsy)
        obs_index_sample_G6_biopsy = unique(obs_index_sample_G6_biopsy)
        
        for(i in 1:length(obs_index_sample_G6_biopsy)){
          future_data = patient_df[1,]
          future_data$visitTimeYears = next_decision_epoch
          future_data$log2psaplus1 = future_outcome_cur_G6$predicted_psa[1, obs_index_sample_G6_biopsy[i]]
          future_data$high_dre = rbinom(1,1, future_outcome_cur_G6$predicted_dre_prob[1, obs_index_sample_G6_biopsy[i]])
          
          future_belief_G6_biopsy = getExpectedFutureOutcomes(fitted_JM, 
                                                            rbind(patient_df, future_data), 
                                                            latest_survival_time = current_decision_epoch,
                                                            earliest_failure_time = earliest_failure_time,
                                                            survival_predict_times =  next_decision_epoch,
                                                            M = N_MCMC_ITER)
          G6_prob = mean(future_belief_G6_biopsy$predicted_surv_prob[1,])
          total_G6 = sum(rbinom(obs_distribution[i], 1, G6_prob))
          
          biopsy_children[[i]] = createNode(rbind(patient_df, future_data), 
                                            next_decision_epoch,
                                            c(G6 = total_G6, G7=obs_distribution[i]-total_G6, AT=0),
                                            current_decision_epoch, current_decision_epoch,
                                            earliest_failure_time, nr_discretized_obs, 
                                            max_decision_epoch)
        }
      }
      
      #For action Biopsy for state is G7
      if(state_distribution[G7] > 0){
        biopsy_children[[length(biopsy_children)]] = list(decision_epoch = next_decision_epoch,
                                                          state_distribution = c(G6=0, G7=0, AT=state_distribution[G7]),
                                                          biopsy_children = NULL, wait_children=NULL)
      }
    }
    
    ####################
    #For wait action
    ####################
    if(state_distribution[G6] > 0){
      if(is.null(future_outcome_cur_G6)){
        future_outcome_cur_G6 = getExpectedFutureOutcomes(fitted_JM, patient_df, 
                                                          latest_survival_time = current_decision_epoch,
                                                          earliest_failure_time = earliest_failure_time,
                                                          dre_predict_times = next_decision_epoch,
                                                          psa_predict_times = next_decision_epoch,
                                                          addRandomError = T, M = N_MCMC_ITER)
        
        #At the moment I resample from the last nr_discretized_obs
        #assuming that nr_discretized_obs represent the distribution of Y well
        obs_index_sample_G6_wait = sample(N_MCMC_ITER:(N_MCMC_ITER-nr_discretized_obs + 1), 
                                          size=state_distribution[G6], replace = T)
        obs_distribution_G6_wait = table(obs_index_sample_G6_wait)
        obs_index_sample_G6_wait = unique(obs_index_sample_G6_wait)
      }
      
    }
    
    if(state_distribution[G7] > 0){  
      future_outcome_G7_wait = getExpectedFutureOutcomes(fitted_JM, patient_df, 
                                                         latest_survival_time = latest_survival_time,
                                                         earliest_failure_time = current_decision_epoch,
                                                         dre_predict_times = next_decision_epoch,
                                                         psa_predict_times = next_decision_epoch,
                                                         addRandomError = T, M = N_MCMC_ITER)
      
      obs_index_sample_G7_wait = sample(N_MCMC_ITER:(N_MCMC_ITER-nr_discretized_obs + 1), 
                                        size=state_distribution[G7], replace = T)
      obs_distribution_G7_wait = table(obs_index_sample_G7_wait)
      obs_index_sample_G7_wait = unique(obs_index_sample_G7_wait)
      
      #We use these later with wait action under G6. These and the ones from G6
      #make the belief at the next level
    }
    
    
  }
  
  return(list(decision_epoch = current_decision_epoch,
              state_distribution = state_distribution,
              biopsy_children = biopsy_children, wait_children = wait_children))
}

selectAction = function(patient_df, current_decision_epoch,
                        latest_survival_time, earliest_failure_time,
                        nr_discretized_obs, max_decision_epoch) {
  
  initial_belief = getExpectedFutureOutcomes(fitted_JM, patient_df, 
                                             latest_survival_time, earliest_failure_time,
                                             survival_predict_times = current_decision_epoch,
                                             M=N_MCMC_ITER)
  
  G6_prob = mean(initial_belief$predicted_surv_prob[1,])
  total_G6 = sum(rbinom(n_scenarios, 1, G6_prob))
  
  despot_tree = createNode(patient_df, current_decision_epoch, 
                           state_distribution = c(G6 = total_G6, G7 = N_DESPOT_SCENARIOS-total_G6, AT=0),
                           latest_biopsy_time = latest_survival_time,
                           latest_survival_time, earliest_failure_time,
                           nr_discretized_obs, max_decision_epoch)
  
  value_despot = computeValue(despot_tree)
  
  return(value_despot)
}

# computeValue = function(DESPOT_TREE){
#   immediate_reward = mean(sapply(DESPOT_TREE, function(x){
#     generateReward(x$state, action)
#   }))
#   
#   future_value = sapply(DESPOT_TREE, function(x){
#     computeValue(x[[action]])
#   })
#   
#   return(immediate_reward + future_value$value)
# }