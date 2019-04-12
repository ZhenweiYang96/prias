createNewDespotNode = function(decision_epoch){
  return(list(decision_epoch = decision_epoch,
              next_decision_epoch = getNextDecisionEpoch(decision_epoch),
              state_dist = c()))
}

simulateScenario = function(scenario_id, belief_id, patient_df, cat_outcomes, latest_survival_time,
                            earliest_failure_time, max_decision_epoch, no_Y){
  
  if(earliest_failure_time==Inf){#i.e. if the previous state_particle was G7
    scenario_signature = paste(DESPOT_TREE[[belief_id]]$decision_epoch,
                               latest_survival_time, 
                               earliest_failure_time, 
                               cat_outcomes, sep = "-")
    if(is.null(DESPOT_BELIEF_CACHE[[scenario_signature]])){
      DESPOT_BELIEF_CACHE[[scenario_signature]] <<- getExpectedFuture_Cat(object = fitted_JM, patient_data = patient_df, 
                                                                          lower_upper_psa_limits = LOWER_UPPER_PSA_LIMITS,
                                                                          latest_survival_time = latest_survival_time,
                                                                          earliest_failure_time = earliest_failure_time, 
                                                                          survival_predict_times = DESPOT_TREE[[belief_id]]$decision_epoch,
                                                                          M = N_MCMC_ITER)$predicted_surv_prob  
    }
    
    p_nocancer = DESPOT_BELIEF_CACHE[[scenario_signature]]
    
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
        if(no_Y == FALSE){
          fut_scenario_signature = paste(DESPOT_TREE[[belief_id]]$next_decision_epoch,
                                         fut_latest_survival_time, 
                                         fut_earliest_failure_time, 
                                         cat_outcomes, sep = "-")
          
          if(is.null(DESPOT_Y_CACHE[[fut_scenario_signature]])){
            DESPOT_Y_CACHE[[fut_scenario_signature]] <<- getExpectedFuture_Cat(object = fitted_JM, patient_data = patient_df, 
                                                                               lower_upper_psa_limits = LOWER_UPPER_PSA_LIMITS,
                                                                               latest_survival_time = fut_latest_survival_time,
                                                                               earliest_failure_time = fut_earliest_failure_time, 
                                                                               Y_predict_times = DESPOT_TREE[[belief_id]]$next_decision_epoch,
                                                                               M = N_MCMC_ITER)$predicted_Y_prob
          }
          
          fut_y_cat_prob = DESPOT_Y_CACHE[[fut_scenario_signature]]
          
          fut_y_cat_id = sample(x = OUTCOME_CAT_NAMES, size = 1, prob = fut_y_cat_prob)
          
          new_row = patient_df[1,]
          new_row$visitTimeYears = DESPOT_TREE[[belief_id]]$next_decision_epoch
          new_row$psa_cat_data = T
          new_row[, c("log2psaplus1", "high_dre")] = OUTCOME_PSA_DRE_CAT[fut_y_cat_id,]
        }else{
          new_row = patient_df[0,]
          fut_y_cat_id = OUTCOME_DUMMY
        }
        
        
        new_belief_id = paste0(belief_id,"; ", action,"-", fut_y_cat_id)
        if(is.null(DESPOT_TREE[[new_belief_id]])){
          DESPOT_TREE[[new_belief_id]] <<- createNewDespotNode(decision_epoch = DESPOT_TREE[[belief_id]]$next_decision_epoch)
        }
        
        simulateScenario(scenario_id, new_belief_id, 
                         rbind(patient_df, new_row), 
                         cat_outcomes = paste(cat_outcomes, fut_y_cat_id, sep=":"),
                         fut_latest_survival_time, fut_earliest_failure_time,
                         max_decision_epoch, no_Y)
      }
    }
  }
}

getOptimalAction = function(belief_id, latest_survival_time, earliest_failure_time, max_decision_epoch,
                            cur_biopsies, max_biopsies, no_Y){
  
  current_decision_epoch = DESPOT_TREE[[belief_id]]$decision_epoch
  
  available_actions = if((current_decision_epoch - latest_survival_time) >= MIN_BIOPSY_GAP & cur_biopsies<max_biopsies){
    ACTION_VECTOR
  }else{
    WAIT
  }
  
  state_dist = DESPOT_TREE[[belief_id]]$state_dist
  total_particles = length(state_dist)
  G6_prob = sum(state_dist==G6) / total_particles
  
  optimal_action = NA
  optimal_reward = -Inf
  for(current_action in available_actions){
    current_reward = G6_prob * getReward(G6, current_action, 
                                         current_decision_epoch, latest_survival_time, 
                                         cur_biopsies) +
      (1-G6_prob) * getReward(G7, current_action, 
                              current_decision_epoch, latest_survival_time, 
                              cur_biopsies)
    
    if(no_Y == FALSE){
      fut_optimal_action_chain = vector("list", length(OUTCOME_CAT_NAMES))
    }else{
      fut_optimal_action_chain = vector("list", 1)
    }
    
    #Exploring the future now
    if(DESPOT_TREE[[belief_id]]$next_decision_epoch <= max_decision_epoch){
      future_cur_biopsies = ifelse(current_action==BIOPSY, cur_biopsies + 1, cur_biopsies)
      future_latest_survival_time = ifelse(current_action==BIOPSY, 
                                           current_decision_epoch, 
                                           latest_survival_time)
      
      #Now it is a bit hard to program for the case with Y and without Y so I am using an inelegant trick
      for(outcome_cat_name in OUTCOME_CAT_NAMES){
        
        if(no_Y == T){
          outcome_cat_name = OUTCOME_DUMMY
        }
        
        future_belief_id = paste0(belief_id,"; ", current_action, "-", outcome_cat_name)
        
        if(!is.null(DESPOT_TREE[[future_belief_id]])){
          future_action_reward = getOptimalAction(future_belief_id, 
                                                  future_latest_survival_time, earliest_failure_time,
                                                  max_decision_epoch, future_cur_biopsies, max_biopsies, no_Y)
          
          fut_optimal_action_chain[[outcome_cat_name]] = future_action_reward$optimal_action_chain
          
          fut_y_prob_empirical = length(DESPOT_TREE[[future_belief_id]]$state_dist) / total_particles
          
          current_reward = current_reward + DISCOUNT_FACTORS[as.character(DESPOT_TREE[[belief_id]]$next_decision_epoch)] * 
            fut_y_prob_empirical * future_action_reward$optimal_reward
        }
        
        #This is a bit inelegent
        if(no_Y == TRUE){
          break
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
  
  return(list(optimal_action=optimal_action, optimal_action_chain=optimal_action_chain,
              optimal_reward=optimal_reward))
}

selectAction = function(patient_df, current_decision_epoch,
                        latest_survival_time, earliest_failure_time,
                        max_decision_epoch, cur_biopsies=0, max_biopsies=Inf,
                        no_Y=F) {
  DESPOT_TREE <<- list()
  DESPOT_BELIEF_CACHE <<- list()
  DESPOT_Y_CACHE <<- list()
  belief_id = "|"
  DESPOT_TREE[[belief_id]] <<- createNewDespotNode(decision_epoch = current_decision_epoch)
  
  for(scenario_id in as.character(1:N_DESPOT_SCENARIOS)){
    simulateScenario(scenario_id, belief_id, patient_df, cat_outcomes = "",
                     latest_survival_time, earliest_failure_time,
                     max_decision_epoch, no_Y)
  }
  
  #Now we traverse through the DESPOT to find the optimal action
  getOptimalAction(belief_id, latest_survival_time, earliest_failure_time, 
                   max_decision_epoch, cur_biopsies, max_biopsies, no_Y)
}