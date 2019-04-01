validateDespotBelief = function(patient_df, current_decision_epoch,
                                latest_survival_time, earliest_failure_time,
                                max_decision_epoch) {
  
  empirical_p_nocancer = rep(NA, length(DESPOT_TREE))
  true_p_nocancer = rep(NA, length(DESPOT_TREE))
  
  biopsy_test_times = BIOPSY_TEST_TIMES[BIOPSY_TEST_TIMES >= current_decision_epoch & 
                                          BIOPSY_TEST_TIMES < max_decision_epoch]
  
  state_dist_length = rep(NA, length(DESPOT_TREE))
  for(i in 1:length(DESPOT_TREE)){
    belief_name = names(DESPOT_TREE)[i]
    decision_epoch = DESPOT_TREE[[i]]$decision_epoch
    state_dist = DESPOT_TREE[[i]]$state_dist
    state_dist_length[i] = length(state_dist)
    empirical_p_nocancer[i] = sum(state_dist==G6) / length(state_dist)
    
    splits = strsplit(belief_name, split = "; ")[[1]]
    
    cur_latest_survival_time = latest_survival_time
    cur_pat_data = patient_df
    
    if(length(splits)>1){
      act_obs = do.call('rbind', strsplit(splits[-1], split="-"))
      
      if(tail(act_obs[,2],1) == OUTCOME_TREATMENT){
        true_p_nocancer[i] = 0
        next
      }
      
      cur_latest_survival_time = max(cur_latest_survival_time,
                                     biopsy_test_times[max(which(act_obs[,1]==BIOPSY))],
                                     na.rm = T)
      
      obs_times = BIOPSY_TEST_TIMES[BIOPSY_TEST_TIMES > current_decision_epoch & 
                                      BIOPSY_TEST_TIMES <= decision_epoch]
      new_obs = OUTCOME_PSA_DRE_CAT[act_obs[,2],]
      cur_pat_data = patient_df[rep(1,nrow(act_obs)),]
      cur_pat_data$psa_cat_data = T
      cur_pat_data$visitTimeYears = obs_times
      
      cur_pat_data[,c("log2psaplus1", "high_dre")] = new_obs
      cur_pat_data = rbind(patient_df, cur_pat_data)
    }
    
    true_p_nocancer[i] = getExpectedFuture_Cat(fitted_JM, cur_pat_data, LOWER_UPPER_PSA_LIMITS,
                                               cur_latest_survival_time,earliest_failure_time,
                                               decision_epoch, M = N_MCMC_ITER)$predicted_surv_prob
  }
  
  return(list(empirical_p_nocancer, true_p_nocancer, state_dist_length))
}

#validation_res = validateDespotBelief(pat_data, 1, 0, Inf, 2)

#DESPOT Value Validation
validateDespotValue = function(belief_id, action_chain, latest_survival_time, 
                               earliest_failure_time, max_decision_epoch,
                               cur_biopsies, max_biopsies){
  
  current_decision_epoch = DESPOT_TREE[[belief_id]]$decision_epoch
  
  current_action = action_chain$Now
  
  state_dist = DESPOT_TREE[[belief_id]]$state_dist
  total_particles = length(state_dist)
  G6_prob = sum(state_dist==G6) / total_particles
  
  optimal_action = NA
  optimal_reward = -Inf
  current_reward = G6_prob * getReward(G6, current_action, 
                                       current_decision_epoch, latest_survival_time, 
                                       cur_biopsies) +
    (1-G6_prob) * getReward(G7, current_action, 
                            current_decision_epoch, latest_survival_time, 
                            cur_biopsies)
  
  #Exploring the future now
  if(DESPOT_TREE[[belief_id]]$next_decision_epoch <= max_decision_epoch){
    future_cur_biopsies = ifelse(current_action==BIOPSY, cur_biopsies + 1, cur_biopsies)
    future_latest_survival_time = ifelse(current_action==BIOPSY, 
                                         current_decision_epoch, 
                                         latest_survival_time)
    
    for(outcome_cat_name in OUTCOME_CAT_NAMES){
      future_belief_id = paste0(belief_id,"; ", current_action, "-", outcome_cat_name)
      
      if(!is.null(DESPOT_TREE[[future_belief_id]])){
        future_action_reward = validateDespotValue(future_belief_id, action_chain$Future[[outcome_cat_name]],
                                                   future_latest_survival_time, earliest_failure_time,
                                                   max_decision_epoch, future_cur_biopsies, max_biopsies)
        
        fut_y_prob_empirical = length(DESPOT_TREE[[future_belief_id]]$state_dist) / total_particles
        
        current_reward = current_reward + DISCOUNT_FACTORS[as.character(DESPOT_TREE[[belief_id]]$next_decision_epoch)] * 
          fut_y_prob_empirical * future_action_reward$optimal_reward
      }
    }
  }
  
  if(current_reward > optimal_reward){
    optimal_action = current_action
    optimal_reward = current_reward
  }
  
  return(list(optimal_action=optimal_action,
              optimal_reward=optimal_reward))
}
