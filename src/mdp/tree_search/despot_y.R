createNewDespotNode = function(decision_epoch){
  return(list(decision_epoch = decision_epoch,
              next_decision_epoch = getNextDecisionEpoch(decision_epoch),
              state_dist = c()))
}

simulateScenario = function(scenario_id, belief_id, patient_df, cat_outcomes, latest_survival_time,
                            earliest_failure_time, max_decision_epoch){
  
  if(earliest_failure_time==Inf){#i.e. if the previous state_particle was G7
    belief_signature = paste(DESPOT_TREE[[belief_id]]$decision_epoch,
                             latest_survival_time, 
                             earliest_failure_time, 
                             cat_outcomes, sep = "-")
    if(is.null(DESPOT_BELIEF_CACHE[[belief_signature]])){
      DESPOT_BELIEF_CACHE[[belief_signature]] <<- getExpectedFuture_Cat(object = fitted_JM, patient_data = patient_df, 
                                                                        lower_upper_psa_limits = LOWER_UPPER_PSA_LIMITS,
                                                                        latest_survival_time = latest_survival_time,
                                                                        earliest_failure_time = earliest_failure_time, 
                                                                        survival_predict_times = DESPOT_TREE[[belief_id]]$decision_epoch,
                                                                        M = N_MCMC_ITER)$predicted_surv_prob  
    }
    
    p_nocancer = DESPOT_BELIEF_CACHE[[belief_signature]]
    
  }else{
    p_nocancer = 0
  }
  
  new_state_particle = sample(c(G6, G7), size = 1, prob = c(p_nocancer, 1-p_nocancer))
  DESPOT_TREE[[belief_id]]$state_dist[scenario_id] <<- new_state_particle
  
  #If future needs to simulated
  if(DESPOT_TREE[[belief_id]]$next_decision_epoch <= max_decision_epoch){
    
    fut_latest_survival_time = ifelse(new_state_particle==G6, 
                                      yes = DESPOT_TREE[[belief_id]]$decision_epoch,
                                      no = latest_survival_time)
    
    fut_earliest_failure_time = ifelse(new_state_particle==G7,
                                       yes = min(earliest_failure_time, DESPOT_TREE[[belief_id]]$decision_epoch),
                                       no = earliest_failure_time)
    
    for(action in ACTION_VECTOR){
      if(new_state_particle==G7 & action==BIOPSY){
        new_belief_id = paste0(belief_id,"; ", action,"-", OUTCOME_TREATMENT)
        
        if(is.null(DESPOT_TREE[[new_belief_id]])){
          DESPOT_TREE[[new_belief_id]] <<- createNewDespotNode(decision_epoch = DESPOT_TREE[[belief_id]]$next_decision_epoch)
        }
        
        #Reached the terminal stage
        DESPOT_TREE[[new_belief_id]]$state_dist[scenario_id] <<- AT
      }else{
        fut_belief_signature = paste(DESPOT_TREE[[belief_id]]$next_decision_epoch,
                                     fut_latest_survival_time, 
                                     fut_earliest_failure_time, 
                                     cat_outcomes, sep = "-")
        
        if(is.null(DESPOT_Y_CACHE[[fut_belief_signature]])){
          DESPOT_Y_CACHE[[fut_belief_signature]] <<- getExpectedFuture_Cat(object = fitted_JM, patient_data = patient_df, 
                                             lower_upper_psa_limits = LOWER_UPPER_PSA_LIMITS,
                                             latest_survival_time = fut_latest_survival_time,
                                             earliest_failure_time = fut_earliest_failure_time, 
                                             Y_predict_times = DESPOT_TREE[[belief_id]]$next_decision_epoch,
                                             M = N_MCMC_ITER)$predicted_Y_prob
        }
        
        fut_y_cat_prob = DESPOT_Y_CACHE[[fut_belief_signature]]
        
        fut_y_cat_id = sample(x = OUTCOME_CAT_NAMES, size = 1, prob = fut_y_cat_prob)
        
        new_row = patient_df[1,]
        new_row$visitTimeYears = DESPOT_TREE[[belief_id]]$next_decision_epoch
        new_row$psa_cat_data = T
        new_row[, c("log2psaplus1", "high_dre")] = OUTCOME_PSA_DRE_CAT[fut_y_cat_id,]
        
        new_belief_id = paste0(belief_id,"; ", action,"-", fut_y_cat_id)
        if(is.null(DESPOT_TREE[[new_belief_id]])){
          DESPOT_TREE[[new_belief_id]] <<- createNewDespotNode(decision_epoch = DESPOT_TREE[[belief_id]]$next_decision_epoch)
        }
        
        simulateScenario(scenario_id, new_belief_id, 
                         rbind(patient_df, new_row), 
                         cat_outcomes = paste(cat_outcomes, fut_y_cat_id, sep=":"),
                         fut_latest_survival_time, fut_earliest_failure_time,
                         max_decision_epoch)
      }
    }
  }
}

selectAction = function(patient_df, current_decision_epoch,
                        latest_survival_time, earliest_failure_time,
                        max_decision_epoch) {
  DESPOT_TREE <<- list()
  DESPOT_BELIEF_CACHE <<- list()
  DESPOT_Y_CACHE <<- list()
  belief_id = "|"
  DESPOT_TREE[[belief_id]] <<- createNewDespotNode(decision_epoch = current_decision_epoch)
  
  for(scenario_id in as.character(1:N_DESPOT_SCENARIOS)){
    simulateScenario(scenario_id, belief_id, patient_df, cat_outcomes = "",
                     latest_survival_time, earliest_failure_time,
                     max_decision_epoch)
  }
}