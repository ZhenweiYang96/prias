generateFuture = function(state, action, patient_df, current_decision_epoch,
                          latest_survival_time, earliest_failure_time){
  next_decision_epoch = getNextDecisionEpoch(current_decision_epoch, action)
  
  if(current_state == G7 & action==BIOPSY){
    future_obs = patient_df[0,]
  }else{
    future_obs = generateObservation(current_decision_epoch, next_decision_epoch,
                                     patient_df, latest_survival_time, 
                                     earliest_failure_time)
  }
  patient_df = rbind(patient_df, future_obs)
  
  if(current_state == G7){
    future_state = ifelse(action==BIOPSY, AT, G7)
  }else{
    future_belief = getExpectedFutureOutcomes(fitted_JM, patient_df,
                                              latest_survival_time, earliest_failure_time,
                                              survival_predict_times = next_decision_epoch,
                                              M=N_MCMC_ITER)
    
    G6_prob = mean(future_belief$predicted_surv_prob[1,])
    
    future_state = sample(x=STATE_VECTOR, size=1, prob=c(G6_prob, 1-G6_prob, 0))
  }
  
  return(list(state = future_state, patient_df = patient_df,
              decision_epoch = next_decision_epoch))
}

createParticle = function(state, patient_df, 
                          current_decision_epoch, 
                          latest_survival_time, earliest_failure_time,
                          max_decision_epoch){
  
  particle = list(state = state, time=current_decision_epoch)
  
  if(state != AT & current_decision_epoch < max_decision_epoch){
    if(current_state==G6) {
      latest_survival_time = current_decision_epoch
    }else if(is.infinite(earliest_failure_time)){
      earliest_failure_time = current_decision_epoch
    }
    
    for(action in ACTION_VECTOR){
      future = generateFuture(state, action, patient_df, current_decision_epoch,
                              latest_survival_time, earliest_failure_time)
      particle[[action]] = createParticle(future$state, future$patient_df,
                                          future$decision_epoch, 
                                          latest_survival_time, earliest_failure_time,
                                          max_decision_epoch)
    }
  }
  
  return(particle)
}

computeValue = function(DESPOT_TREE){
  immediate_reward = mean(sapply(DESPOT_TREE, function(x){
    generateReward(x$state, action)
  }))
  
  future_value = sapply(DESPOT_TREE, function(x){
    computeValue(x[[action]])
  })
  
  return(immediate_reward + future_value$value)
}

selectAction = function(patient_df, n_scenarios,
                        current_decision_epoch,
                        latest_survival_time, earliest_failure_time,
                        discount_factors, max_decision_epoch) {
  
  initial_belief = getExpectedFutureOutcomes(fitted_JM, patient_df,
                                             latest_survival_time, earliest_failure_time,
                                             survival_predict_times = current_decision_epoch,
                                             M=N_MCMC_ITER)
  
  G6_prob = mean(initial_belief$predicted_surv_prob[1,])
  
  #First we build the DESPOT
  sampled_states = sample(x=STATE_VECTOR, size=n_scenarios, replace = T, 
                          prob=c(G6_prob, 1-G6_prob, 0))
  
  #Assuming DESPOT_TREE is an empty list
  for(i in 1:n_scenarios){
    DESPOT_TREE[[i]] <<- createParticle(sampled_states[i], patient_df, 
                                        current_decision_epoch,
                                        latest_survival_time, earliest_failure_time,
                                        max_decision_epoch)
  }
  
  value_despot = computeValue(DESPOT_TREE)
  
  return(value_despot)
}