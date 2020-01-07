surv_formula = ~ age2 + gender + BElength_cat + esophagitis
long_formula = ~ 1 + time + age2 + gender + BElength_cat + esophagitis

gammas = jmfit_barrett_value_expit$statistics$postMeans$gammas
beta_sox2 = jmfit_barrett_value_expit$statistics$postMeans$betas1
beta_p53 = jmfit_barrett_value_expit$statistics$postMeans$betas2
beta_LGD = jmfit_barrett_value_expit$statistics$postMeans$betas3

generateBaselineData = function(n_sub){
  #Generate gender
  gender = relevel(factor(sample(c("male", "female"), size = n_sub,replace = T, c(0.7337559, 1-0.7337559))), 
                   ref = "male")
  maleMeanAge2 = mean(data.id$age2[data.id$gender=="male"])
  maleSDAge2 = sd(data.id$age2[data.id$gender=="male"])
  femaleMeanAge2 = mean(data.id$age2[data.id$gender=="female"])
  femaleSDAge2 = sd(data.id$age2[data.id$gender=="female"])
  
  age2 = sapply(gender, FUN = function(x){
    if(x=="male"){
      return(rnorm(1, maleMeanAge2, maleSDAge2))
    }else{
      return(rnorm(1, femaleMeanAge2, femaleSDAge2))
    }
  })
  id = 1:n_sub
  
  BElength_cat_prob = plogis(3.2072 + rnorm(n_sub, mean=0, sd=3.992))
  esophagitis_log_odds_notime = -3.69772 + rnorm(n_sub, mean=0, sd = 1.829)
  
  #every week
  #exo_times = seq(0, MAX_FAIL_TIME, length.out=MAX_FAIL_TIME * 50)
  exo_times = 0
  patient_data = data.frame(id=rep(id, each=length(exo_times)), 
                            age2 = rep(age2, each=length(exo_times)),
                            gender = rep(gender, each=length(exo_times)),
                            time = rep(exo_times, n_sub))
  patient_data$BElength_cat = sapply(rep(BElength_cat_prob, each=length(exo_times)), rbinom, n=1, size=1)
  patient_data$esophagitis = sapply(plogis(rep(esophagitis_log_odds_notime, each=length(exo_times)) - 0.09237 * patient_data$time), 
                                    rbinom, n=1, size=1)
  ###########################
  # Generate the longitudinal dataset
  ###########################
  D = jmfit_barrett_value_expit$statistics$postMeans$D
  b <- mvrnorm(n_sub, mu = rep(0, nrow(D)), D)
  
  gammas = jmfit_barrett_value_expit$statistics$postMeans$gammas
  beta_sox2 = jmfit_barrett_value_expit$statistics$postMeans$betas1
  beta_p53 = jmfit_barrett_value_expit$statistics$postMeans$betas2
  beta_LGD = jmfit_barrett_value_expit$statistics$postMeans$betas3
  
  hazardFunc = function (time, patient_id) {
    pat_data = patient_data[rep(patient_id, length(time)),]
    pat_data$time = time
    
    wGamma = model.matrix(surv_formula, pat_data)[,-1] %*% gammas
    baselinehazard = exp(splineDesign(jmfit_barrett_value_expit$control$knots, time, 
                                      ord = jmfit_barrett_value_expit$control$ordSpline, outer.ok = T) %*% jmfit_barrett_value_expit$statistics$postMeans$Bs_gammas)
    
    long_X = model.matrix(long_formula, pat_data)
    b_subject = b[patient_id,]
    true_sox2_prob = plogis(long_X %*% beta_sox2 + b_subject[1])
    true_p53_prob = plogis(long_X %*% beta_p53 + b_subject[2])
    true_LGD_prob = plogis(long_X %*% beta_LGD + b_subject[3])
    
    y_Alpha = cbind(true_sox2_prob, true_p53_prob, true_LGD_prob) %*% jmfit_barrett_value_expit$statistics$postMeans$alphas
    
    baselinehazard * exp(wGamma + y_Alpha)
  }
  
  invSurvival <- function (t, u, patient_id) {
    #gk_weights = 0.5 * t * JMbayes:::gaussKronrod()$wk
    #gk_points = 0.5 * t * JMbayes:::gaussKronrod()$sk + 0.5 * t
    log(u) + integrate(hazardFunc, lower = 0.001, upper = t, patient_id, stop.on.error = T)$value
    #log(u) + sum(gk_weights * hazardFunc(gk_points, patient_id))
  }
  
  survivalFunc <- function (t, patient_id) {
    #gk_weights = 0.5 * t * JMbayes:::gaussKronrod()$wk
    #gk_points = 0.5 * t * JMbayes:::gaussKronrod()$sk + 0.5 * t
    exp(-integrate(hazardFunc, lower = 0.001, upper = t, patient_id, stop.on.error = T)$value)
    #exp(-sum(gk_weights * hazardFunc(gk_points, patient_id)))
  }
  
  fastSurvivalFunc <- function (t, patient_id) {
    gk_weights = 0.5 * t * JMbayes:::gaussKronrod()$wk
    gk_points = 0.5 * t * JMbayes:::gaussKronrod()$sk + 0.5 * t
    #exp(-integrate(hazardFunc, lower = 0.001, upper = t, patient_id, stop.on.error = F)$value)
    exp(-sum(gk_weights * hazardFunc(gk_points, patient_id)))
  }
  
  pSurvTime = function(surv_prob, patient_id){
    Low = 1e-05
    Up <- MAX_FAIL_TIME
    
    #print(patient_id)
    if(fastSurvivalFunc(MAX_FAIL_TIME, patient_id) > surv_prob){
      return(MAX_FAIL_TIME)
    }
    
    Root <- try(uniroot(invSurvival, interval = c(Low, Up), 
                        u = surv_prob, patient_id=patient_id, tol = 1e-3)$root, TRUE)
    if(inherits(Root, "try-error")){
      print(Root)
      return(MAX_FAIL_TIME)
    }else{
      return(Root)
    }
  }
  
  survProbs = runif(n=n_sub, 0, 1)
  ct = makeCluster(max_cores)
  registerDoParallel(ct)
  prog_time = foreach(i=1:n_sub, .combine='c', 
                      .export = c('MAX_FAIL_TIME','patient_data',
                                  'surv_formula','long_formula',
                                  'jmfit_barrett_value_expit'),
                      .packages = c("splines", "JMbayes")) %dopar%{
                        pSurvTime(survProbs[i], i)
                      }
  stopCluster(ct)
  
  patient_data$prog_time = rep(prog_time, each = length(exo_times))
  patient_data$status = ifelse(patient_data$prog_time==MAX_FAIL_TIME, 0, 1)
  print(paste("Event rate:", 100 * sum(patient_data$status)/n_sub, "%"))
  
  # long_times = c(0, 3/12, 6/12, 9/12, seq(1, MAX_FAIL_TIME, 1))
  # patient_data_long = patient_data[rep(1:n_sub, each=length(long_times)),]
  # patient_data_long$time = rep(long_times, n_sub)
  # 
  # longX_matrix = model.matrix(long_formula, patient_data_long)
  # 
  # #Now we add longitudinal data to this patient
  # patient_data_long$sox2 = sapply(plogis(longX_matrix %*% beta_sox2 + rep(b[,1], each=length(long_times))),
  #                            rbinom, n=1, size=1)
  # patient_data_long$p53 = sapply(plogis(longX_matrix %*% beta_p53 + rep(b[,2], each=length(long_times))),
  #                           rbinom, n=1, size=1)
  # patient_data_long$LGD = sapply(plogis(longX_matrix %*% beta_LGD + rep(b[,3], each=length(long_times))),
  #                           rbinom, n=1, size=1)
  # 
  # patient_data_long = patient_data_long[patient_data_long$time <= patient_data_long$prog_time,]
  
  retList = list(patient_data = patient_data, b = b)
  return(retList)
}

fitJointModel = function(patient_data, patient_data_long, n_training, 
                         cens_start_time = 0, cens_end_time = 30){
  
  cox_part = coxph(Surv(prog_time, status)~age2 + gender + BElength_cat + esophagitis, 
                   data=patient_data, model = T, x = T)
  
  mvglmer_part <- mvglmer(list(
    sox2 ~ time + age2 + gender + BElength_cat + esophagitis + (1 | id),
    p53  ~ time + age2 + gender + BElength_cat + esophagitis + (1 | id),
    LGD  ~ time + age2 + gender + BElength_cat + esophagitis + (1 | id)
  ), data = patient_data_long, families = list(binomial, binomial, binomial), 
  engine = "JAGS", control=list(n.processors=max_cores))
  
  fitted_jm = mvJointModelBayes(mvglmer_part, cox_part, 
                                timeVar = "time", 
                                transFuns = tFuns, 
                                priors = list(shrink_gammas = TRUE, shrink_alphas = TRUE), 
                                control =  list(n_cores=max_cores))
  
  return(list(patient_data = patient_data, patient_data_long=patient_data_long,
              cox_part = cox_part, mvglmer_part = mvglmer_part,
              fitted_jm = fitted_jm))
}

runFixedSchedule = function(test_data){
  ct = makeCluster(max_cores)
  registerDoParallel(ct)
  nb_offset = foreach(i=unique(test_data$id), .combine='rbind', 
                      .export=c("simulateProtocol", "MAX_FAIL_TIME"),
                      .packages = c("splines", "JMbayes")) %dopar%{
                        test_data.i = simulateProtocol(test_data[test_data$id == i,])
                        
                        test_data.i$time[nrow(test_data.i)] = min(MAX_FAIL_TIME,
                                                                  test_data.i$time[nrow(test_data.i)])
                        
                        offset = tail(test_data.i$time,1) - test_data.i$prog_time[1]
                        return(c("nb"=nrow(test_data.i), "offset"=offset))
                      }
  stopCluster(ct)
  
  return(nb_offset)
}

runRiskBasedSchedule = function(fitted_jm, test_data, threshold, min_gap=0.5){
  ct = makeCluster(max_cores)
  registerDoParallel(ct)
  nb_offset = foreach(i=unique(test_data$id), .combine='rbind', 
                      .export=c("simulateRiskBased", "MAX_FAIL_TIME"),
                      .packages = c("splines", "JMbayes")) %dopar%{
                        test_data.i = simulateRiskBased(fitted_jm, test_data[test_data$id == i,], threshold, min_gap)
                        
                        test_data.i$time[nrow(test_data.i)] = min(MAX_FAIL_TIME,
                                                                  test_data.i$time[nrow(test_data.i)])
                        
                        offset = tail(test_data.i$time,1) - test_data.i$prog_time[1]
                        return(c("nb"=nrow(test_data.i), "offset"=offset))
                      }
  stopCluster(ct)
  
  return(nb_offset)
}

#Send 1 row of data
simulateProtocol = function(patient_data){
  available_time_list = patient_data$time
  
  prog_time = patient_data$prog_time[1]
  
  #Choose 6 months of data
  ret_data = patient_data[which.min(abs(available_time_list-0)),]
  
  latest_gap = 3
  next_time = 0.5
  repeat{
    cur_time = next_time
    new_data = patient_data[which.min(abs(available_time_list-cur_time)),]
    new_data$time = cur_time
    ret_data = rbind(ret_data, new_data)
    
    if(new_data$LGD == 1){
      if(latest_gap==3){
        next_time = cur_time + 0.5
      }else{
        next_time = cur_time + 1
      }
    }else{
      if(latest_gap==0.5){
        next_time = cur_time + 1
      }else{
        next_time = cur_time + 3
      }
    }
    
    latest_gap = next_time - cur_time
    
    if(cur_time > min(prog_time, MAX_FAIL_TIME, na.rm = T)){
      if(!is.na(prog_time)){
        ret_data$LGD[nrow(ret_data)] = NA
      }
      break
    }
  }
  
  return(ret_data)
}

simulateRiskBased = function(fitted_jm, patient_data, threshold, min_gap=0.5){
  pDynSurv = function(survProb, patientDs){
    invDynSurvLastTime <- function (time, u, patientDs, maxRiskDt) {
      u - survfitJM(fitted_jm, patientDs, survTimes = time)$summaries[[1]][1, "Mean"]
    }
    
    Low = max(patientDs$time) + 1e-05
    Up <- MAX_FAIL_TIME
    
    Root <- try(uniroot(invDynSurvLastTime, interval = c(Low, Up), 
                        u = survProb, patientDs = patientDs)$root, TRUE)
    
    if(inherits(Root, "try-error")){
      return(MAX_FAIL_TIME)
    }else{
      return(Root)
    }
  }
  
  ############## main algo starts here ############
  available_time_list = patient_data$time
  prog_time = patient_data$prog_time[1]
  
  #Choose 6 months of data
  ret_data = patient_data[sapply(c(0,0.5), FUN = function(x){
    which.min(abs(available_time_list - x))
  }),]
  
  ret_data$time = c(0,0.5)
  
  while(tail(ret_data$time, 1) < min(prog_time, MAX_FAIL_TIME, na.rm = T)){
    new_time = pDynSurv(survProb = threshold, patientDs = ret_data)
    
    if(new_time - tail(ret_data$time, 1) <= min_gap){
      new_time = tail(ret_data$time, 1) + min_gap
    }
    
    new_data = patient_data[which.min(abs(available_time_list-new_time)),]
    new_data$time = new_time
    ret_data = rbind(ret_data, new_data)
  }
  
  if(!is.na(prog_time)){
    ret_data$LGD[nrow(ret_data)] = NA
  }
  
  return(ret_data)
}
