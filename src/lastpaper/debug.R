seeds = 2019:2023

combined_df = data.frame(seed=numeric(), nb1_to_delay=numeric(),delay_limit=numeric(),
                         MAX_FOLLOW_UP=numeric(), P_ID=numeric(), age=numeric(),
                         progression_time=numeric(), survProbs=numeric(),nb=numeric(), delay=numeric())

for(sim_seed in seeds){
  # for(nb1_to_delay in c(0.25,0.5,1)){
  for(nb1_to_delay in c(1)){
    #for(delay_limit in c(0.5,1, Inf)){
    for(delay_limit in c(0.05,0.1, 0.15, 0.25, 0.5)){
      for(MAX_FOLLOW_UP in c(6,10)){
        tt=try(load(file = paste("C:/Users/838035/threshold_limit/schedule_res", 
                                 sim_seed, nb1_to_delay, delay_limit, MAX_FOLLOW_UP, ".Rdata", 
                                 sep = "_")), silent = T)
        if(!inherits(tt, 'try-error')){
          newdf = data.frame(seed=sim_seed, nb1_to_delay=nb1_to_delay, delay_limit=delay_limit,
                             MAX_FOLLOW_UP=MAX_FOLLOW_UP, P_ID=sim_res$testData$testDs.id$P_ID,
                             age=sim_res$testData$testDs.id$age, progression_time=sim_res$testData$testDs.id$progression_time,
                             survProbs=sim_res$testData$testDs.id$survProbs, nb=sim_res$testData$testDs.id[,13],
                             delay=sim_res$testData$testDs.id[,14])
          combined_df = rbind(combined_df, newdf)
        }
      }
    }
  }
}

personalized_df = combined_df
personalized_df$progression_time = pmin(personalized_df$MAX_FOLLOW_UP, personalized_df$progression_time)
personalized_df$delay[personalized_df$MAX_FOLLOW_UP == personalized_df$progression_time] = 0

#For prias
getPRIASSchedule = function(patient_data, horizon=10){
  
  #making schedule now
  PSA_CHECK_UP_TIME = c(seq(0, 2, 0.25), seq(2.5, 10, 0.5))
  visit_schedule = PSA_CHECK_UP_TIME[PSA_CHECK_UP_TIME >= 1 & PSA_CHECK_UP_TIME <=horizon]
  
  obs_psa = 2^(patient_data$log2psaplus1) - 1
  obs_psa = pmax(obs_psa, 0.01)
  obs_psa_times = patient_data$year_visit
  
  progression_time = min(horizon, patient_data$progression_time)
  
  nb = 0
  delay = -Inf
  fixed_schedule = c(1, 4, 7, 10, 15)
  latest_biopsy_time = 0
  for(i in 1:length(visit_schedule)){
    
    if(latest_biopsy_time >= fixed_schedule[1]){
      fixed_schedule = fixed_schedule[-1]
    }
    
    #in some cases two PSA at current visit time 
    log2psa = log(obs_psa[obs_psa_times <= visit_schedule[i]], base = 2)
    year_visit = obs_psa_times[obs_psa_times <= visit_schedule[i]]
    
    psa_dt = 1/(lm(log2psa~year_visit)$coefficients[2])
    
    if((visit_schedule[i] - latest_biopsy_time) >= 1){
      #if switch to annual schedule
      if(psa_dt>=0 & psa_dt<=10){
        latest_biopsy_time = visit_schedule[i]
        nb = nb + 1
      } else if(visit_schedule[i] >= fixed_schedule[1]){
        latest_biopsy_time = visit_schedule[i]
        nb = nb+1
      }
    }
    
    delay = latest_biopsy_time - progression_time
    if(delay >= 0){
      break
    }
  }
  
  if(delay<0){
    delay = horizon - progression_time
    nb = nb + 1
  }
  return(c(nb, delay))
}

combined_df = data.frame(seed=numeric(), schedule=character(),
                         MAX_FOLLOW_UP=numeric(), P_ID=numeric(), age=numeric(),
                         progression_time=numeric(), survProbs=numeric(),nb=numeric(), delay=numeric())

for(sim_seed in seeds){
  for(MAX_FOLLOW_UP in c(6,10)){
    #for(MAX_FOLLOW_UP in c(10)){
    load(file = paste0("Rdata/lastpaper/sims/sim_seed_",sim_seed,".Rdata"))
    
    yearly = seq(1, MAX_FOLLOW_UP, 1)
    biennial = seq(2, MAX_FOLLOW_UP, 2)
    biennial_after1 = c(1, biennial)
    triennial = unique(c(seq(3, MAX_FOLLOW_UP, 3), MAX_FOLLOW_UP))
    triennial_after1 = c(1, triennial)
    
    nb_delay = t(sapply(sim_res$testData$testDs.id$progression_time, FUN = function(x){
      x=ifelse(x>MAX_FOLLOW_UP, MAX_FOLLOW_UP, no = x)
      nb = which(yearly >= x)[1]
      return(c(nb, yearly[nb]-x))
    }))
    
    newdf = data.frame(seed=sim_seed, schedule="Yearly",
                       MAX_FOLLOW_UP=MAX_FOLLOW_UP,
                       P_ID=sim_res$testData$testDs.id$P_ID,
                       age=sim_res$testData$testDs.id$age,
                       progression_time=sim_res$testData$testDs.id$progression_time,
                       survProbs=sim_res$testData$testDs.id$survProbs,
                       nb=nb_delay[,1],
                       delay=nb_delay[,2])
    
    combined_df = rbind(combined_df, newdf)
    
    nb_delay = t(sapply(sim_res$testData$testDs.id$progression_time, FUN = function(x){
      x=ifelse(x>MAX_FOLLOW_UP, MAX_FOLLOW_UP, no = x)
      nb = which(biennial >= x)[1]
      return(c(nb, biennial[nb]-x))
    }))
    
    newdf = data.frame(seed=sim_seed, schedule="Biennial",
                       MAX_FOLLOW_UP=MAX_FOLLOW_UP,
                       P_ID=sim_res$testData$testDs.id$P_ID,
                       age=sim_res$testData$testDs.id$age,
                       progression_time=sim_res$testData$testDs.id$progression_time,
                       survProbs=sim_res$testData$testDs.id$survProbs,
                       nb=nb_delay[,1],
                       delay=nb_delay[,2])
    
    combined_df = rbind(combined_df, newdf)
    
    nb_delay = t(sapply(sim_res$testData$testDs.id$progression_time, FUN = function(x){
      x=ifelse(x>MAX_FOLLOW_UP, MAX_FOLLOW_UP, no = x)
      nb = which(biennial_after1 >= x)[1]
      return(c(nb, biennial_after1[nb]-x))
    }))
    
    newdf = data.frame(seed=sim_seed, schedule="Biennial_1",
                       MAX_FOLLOW_UP=MAX_FOLLOW_UP,
                       P_ID=sim_res$testData$testDs.id$P_ID,
                       age=sim_res$testData$testDs.id$age,
                       progression_time=sim_res$testData$testDs.id$progression_time,
                       survProbs=sim_res$testData$testDs.id$survProbs,
                       nb=nb_delay[,1],
                       delay=nb_delay[,2])
    
    combined_df = rbind(combined_df, newdf)
    
    nb_delay = t(sapply(sim_res$testData$testDs.id$progression_time, FUN = function(x){
      x=ifelse(x>MAX_FOLLOW_UP, MAX_FOLLOW_UP, no = x)
      nb = which(triennial >= x)[1]
      return(c(nb, triennial[nb]-x))
    }))
    
    newdf = data.frame(seed=sim_seed, schedule="Triennial",
                       MAX_FOLLOW_UP=MAX_FOLLOW_UP,
                       P_ID=sim_res$testData$testDs.id$P_ID,
                       age=sim_res$testData$testDs.id$age,
                       progression_time=sim_res$testData$testDs.id$progression_time,
                       survProbs=sim_res$testData$testDs.id$survProbs,
                       nb=nb_delay[,1],
                       delay=nb_delay[,2])
    
    combined_df = rbind(combined_df, newdf)
    
    nb_delay = t(sapply(sim_res$testData$testDs.id$progression_time, FUN = function(x){
      x=ifelse(x>MAX_FOLLOW_UP, MAX_FOLLOW_UP, no = x)
      nb = which(triennial_after1 >= x)[1]
      return(c(nb, triennial_after1[nb]-x))
    }))
    
    newdf = data.frame(seed=sim_seed, schedule="Triennial_1",
                       MAX_FOLLOW_UP=MAX_FOLLOW_UP,
                       P_ID=sim_res$testData$testDs.id$P_ID,
                       age=sim_res$testData$testDs.id$age,
                       progression_time=sim_res$testData$testDs.id$progression_time,
                       survProbs=sim_res$testData$testDs.id$survProbs,
                       nb=nb_delay[,1],
                       delay=nb_delay[,2])
    
    combined_df = rbind(combined_df, newdf)
    
    nb_delay = by(data = sim_res$testData$testDs, INDICES = sim_res$testData$testDs$P_ID,
                  FUN = getPRIASSchedule, horizon = MAX_FOLLOW_UP)
    
    nb_delay = matrix(unlist(nb_delay), ncol = 2, byrow = T)
    
    newdf = data.frame(seed=sim_seed, schedule="PRIAS",
                       MAX_FOLLOW_UP=MAX_FOLLOW_UP, 
                       P_ID=sim_res$testData$testDs.id$P_ID,
                       age=sim_res$testData$testDs.id$age, 
                       progression_time=sim_res$testData$testDs.id$progression_time,
                       survProbs=sim_res$testData$testDs.id$survProbs, 
                       nb=nb_delay[,1],
                       delay=nb_delay[,2])
    
    combined_df = rbind(combined_df, newdf)
    
  }
}

combined_df$progression_time = pmin(combined_df$MAX_FOLLOW_UP, combined_df$progression_time)
