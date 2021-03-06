fixed_psaFormula = ~ 1 + I(Age - 70) +  I((Age - 70)^2) + ns(visitTimeYears, knots=c(0.1, 0.7, 4), Boundary.knots=c(0, 5.42))
random_psaFormula = ~ 1 + ns(visitTimeYears, knots=c(0.1, 0.7, 4), Boundary.knots=c(0, 5.42))
fixed_random_psaSlopeFormula = ~ 0 + dns(visitTimeYears, knots=c(0.1, 0.7, 4), Boundary.knots=c(0, 5.42))

fixed_dreFormula = ~ 1 + I(Age - 70) +  I((Age - 70)^2) + visitTimeYears
random_dreFormula = ~ 1 + visitTimeYears

survivalFormula = ~ 0 + I(Age - 70) + I((Age - 70)^2)

wk = JMbayes:::gaussKronrod()$wk
sk = JMbayes:::gaussKronrod()$sk

#betas_psa dimension should be p x 1, or p x #MCMC,
#b_psa dimension should be #b x 1, or #b x #MCMC
psaXbetaZb = function(Age, visitTimeYears, betas_psa, b_psa){
  df = data.frame(Age, visitTimeYears)
  model.matrix(fixed_psaFormula, df) %*% betas_psa + 
    model.matrix(random_psaFormula, df) %*% b_psa
}

#betas_psa_slope dimension should be p x 1, or p x #MCMC,
#b_psa_slope dimension should be #b x 1, or #b x #MCMC
psaSlopeXbetaZb = function(Age, visitTimeYears, betas_psa_slope, b_psa_slope){
  X_Z_matrix = model.matrix(fixed_random_psaSlopeFormula, data.frame(Age, visitTimeYears))
  X_Z_matrix %*% betas_psa_slope + X_Z_matrix %*% b_psa_slope
}

#betas_dre dimension should be p x 1, or p x #MCMC,
#b_dre dimension should be #b x 1, or #b x #MCMC
dreLogOddsXbetaZb = function(Age, visitTimeYears, betas_dre, b_dre){
  df = data.frame(Age, visitTimeYears)
  model.matrix(fixed_dreFormula, df) %*% betas_dre +
    model.matrix(random_dreFormula, df) %*% b_dre
}

cum_hazard = function(b, obs_exp_Wgamma, GQ_hazard_data, alphas){
  GQ_hazard_true_dreLogOdds = GQ_hazard_data$GQ_dre_Xbeta + GQ_hazard_data$GQ_dre_Z %*% b[1:2]  
  GQ_hazard_true_psa = GQ_hazard_data$GQ_psa_value_Xbeta + GQ_hazard_data$GQ_psa_value_Z %*% b[3:7]
  GQ_hazard_true_psaSlope = GQ_hazard_data$GQ_psa_slope_Xbeta + GQ_hazard_data$GQ_psa_slope_Z %*% b[4:7]
  
  GQ_hazard_longitudinal_part = exp(GQ_hazard_true_dreLogOdds * alphas[1] + GQ_hazard_true_psa * alphas[2] + GQ_hazard_true_psaSlope * alphas[3])
  
  cumulativeHazard = sum(GQ_hazard_data$GQ_weights_baselinehazard * obs_exp_Wgamma * GQ_hazard_longitudinal_part)
  return(cumulativeHazard)
}

log_p_T = function(b, latest_survival_time, earliest_failure_time, obs_exp_Wgamma, GQ_hazard_fixed_part, alphas){
  cum_hazard_latest_survival_time = if(latest_survival_time==0){
    0
  }else{
    cum_hazard(b, obs_exp_Wgamma, GQ_hazard_fixed_part$latest_survival_time_data, alphas)
  }
  
  cum_hazard_earliest_failure_time = if(earliest_failure_time==Inf){
    Inf
  }else{
    cum_hazard(b, obs_exp_Wgamma, GQ_hazard_fixed_part$earliest_failure_time_data, alphas)
  }
  
  return(log(exp(-cum_hazard_latest_survival_time) - exp(-cum_hazard_earliest_failure_time)))
}

#Just make sure randEff is a vector not a matrix
log_numerator_bayesrule = function(b, data){
  
  #Prior contribution
  randomEff_prior_Contrib = -0.5 * (b %*% data$inv_D %*% b)
  
  #Longitudinal DRE contribution
  fitted_probHighDre = plogis(data$Xbeta_Wgamma_h0BsGamma$obsDre_Xbeta + data$obsDre_Z %*% b[1:2])
  dre_contrib = sum(data$obs_high_dre * log(fitted_probHighDre) + 
                      (1 - data$obs_high_dre) * log(1 - fitted_probHighDre))
  
  #Longitudinal PSA contribution
  Y_fitted_psa = c(data$Xbeta_Wgamma_h0BsGamma$obsPsa_Xbeta + data$obsPsa_Z %*% b[3:7])
  
  Y_psa_meanShift=data$obs_log2psaplus1 - Y_fitted_psa
  meanShift_sigmaMatrix_meanShift = Y_psa_meanShift %*% data$invSigmaMatrix %*% Y_psa_meanShift
  
  psa_contrib = if(data$psaDist=="normal"){
    -0.5 * meanShift_sigmaMatrix_meanShift
  }else{
    #T distribution
    -0.5 * (data$df + length(Y_psa_meanShift)) * log(1 + meanShift_sigmaMatrix_meanShift/data$df)
  }
  
  lower = data$lower_psa_limits[data$obs_log2psaplus1_category]
  upper = data$upper_psa_limits[data$obs_log2psaplus1_category]
  Y_fitted_psa_forcat = c(data$Xbeta_Wgamma_h0BsGamma$obsPsa_cat_Xbeta + data$obsPsa_cat_Z %*% b[3:7])
  
  psa_category_contrib = if(data$psaDist=="normal"){
    sum(sapply(1:length(data$obs_log2psaplus1_category), FUN = function(psa_cat_data_row){
      psa_cat = data$obs_log2psaplus1_category[psa_cat_data_row]
      visitTimeYears = data$obs_psa_cat_visitTimeYears[psa_cat_data_row] 
      psa_cat_quad = data$psa_cat_quad[[as.character(visitTimeYears)]][[psa_cat]]
      fitted_psa_forcat = Y_fitted_psa_forcat[psa_cat_data_row]
      
      if(psa_cat==1){
        #psa_cat_quad[2] is the upper limit, lower limit is -Inf
        pnorm(psa_cat_quad[2], mean=fitted_psa_forcat, 
              sd=data$sigma, log.p = T, lower.tail = T)
      }else if(psa_cat == data$max_psa_cat){
        #psa_cat_quad[1] is the lower limit, upper limit is Inf
        pnorm(psa_cat_quad[1], mean=fitted_psa_forcat, 
              sd=data$sigma, log.p = T, lower.tail = F)
      }else{
        log(sum(psa_cat_quad$weights * 
                  sapply(psa_cat_quad$points, dnorm, 
                         mean=fitted_psa_forcat, sd=data$sigma)))
      }
    }))
  }
  #Time to event contribution
  time_to_event_contrib = log_p_T(b, data$latest_survival_time, data$earliest_failure_time, data$Xbeta_Wgamma_h0BsGamma$obs_exp_Wgamma,
                                  data$Xbeta_Wgamma_h0BsGamma$GQ_hazard_fixed_part, data$alphas)
  
  return(randomEff_prior_Contrib + dre_contrib + psa_contrib + psa_category_contrib + time_to_event_contrib)
}

getGaussianQuadWeightsPoints = function(lower_upper_limit){
  p1 = 0.5 * sum(lower_upper_limit)
  p2 = 0.5 * diff(lower_upper_limit)
  
  return(list(weights = p2*wk, points = p2 * sk + p1))
}

#When M=0 it becomes empirical bayes
get_b_fullBayes = function(object, patient_data, lower_upper_psa_limits,
                           latest_survival_time, earliest_failure_time=Inf,
                           scale = 1.6, psaDist = "normal", TdistDf=3, M=200,
                           optim_methods = c("BFGS", "L-BFGS-B", "CG")){
  
  if(is.null(patient_data$psa_cat_data)){
    patient_data$psa_cat_data = F
  }
  
  #Obs data X and Z matrices for longitudinal part
  patient_data_dre = patient_data[!is.na(patient_data$high_dre), c("Age",  "visitTimeYears", "high_dre")]
  
  
  patient_data_psa = patient_data[!is.na(patient_data$log2psaplus1), c("Age",  "visitTimeYears", "log2psaplus1", "psa_cat_data")]
  
  patient_data_psa_cat = patient_data_psa[patient_data_psa$psa_cat_data==T,]
  patient_data_psa = patient_data_psa[patient_data_psa$psa_cat_data==F,]
  
  if(nrow(patient_data_psa_cat)==0){
    patient_data_psa_cat = patient_data_psa[1,]
    patient_data_psa_cat$log2psaplus1 = 1
    patient_data_psa_cat$visitTimeYears = 0
    lower_upper_psa_limits = list("0"=list(c(-Inf, Inf)))
  }
  
  patient_age = patient_data$Age[1]
  
  obsPsa_X = model.matrix(fixed_psaFormula, patient_data_psa)
  obsPsa_cat_X = model.matrix(fixed_psaFormula, patient_data_psa_cat)
  obsPsa_Z = model.matrix(random_psaFormula, patient_data_psa)
  obsPsa_cat_Z = model.matrix(random_psaFormula, patient_data_psa_cat)
  obsDre_X = model.matrix(fixed_dreFormula, patient_data_dre)
  obsDre_Z = model.matrix(random_dreFormula, patient_data_dre)
  
  obs_W = model.matrix(survivalFormula, data = data.frame(Age = patient_age))
  
  createSurvGQMatrices = function(lower_upper_limit, patient_age, baseline_hazard_knots, baseline_hazard_ordSpline){
    if(sum(lower_upper_limit) != Inf & diff(lower_upper_limit)!=0){
      wp = getGaussianQuadWeightsPoints(lower_upper_limit)
      GQ_h0_matrix =  splineDesign(knots = baseline_hazard_knots, ord = baseline_hazard_ordSpline,
                                   x = wp$points, outer.ok = T)
      
      GQ_patient_df = data.frame(Age = patient_age, visitTimeYears = wp$points)
      
      GQ_psa_value_X = model.matrix(fixed_psaFormula, GQ_patient_df)
      GQ_psa_value_Z = model.matrix(random_psaFormula, GQ_patient_df)
      GQ_psa_slope_X_Z = model.matrix(fixed_random_psaSlopeFormula, GQ_patient_df)
      GQ_dre_X = model.matrix(fixed_dreFormula, GQ_patient_df)
      GQ_dre_Z = model.matrix(random_dreFormula, GQ_patient_df)
      
      return(list(GQ_weights = wp$weights,
                  GQ_h0_matrix = GQ_h0_matrix,
                  GQ_psa_value_X = GQ_psa_value_X,
                  GQ_psa_value_Z = GQ_psa_value_Z,
                  GQ_psa_slope_X_Z = GQ_psa_slope_X_Z,
                  GQ_dre_X = GQ_dre_X,
                  GQ_dre_Z = GQ_dre_Z))
    }else{
      return(NULL)
    }
  }
  
  GQ_hazard_X_Z_W_h0 = lapply(list("latest_survival_time_data"=c(0,latest_survival_time),
                                   "earliest_failure_time_data"=c(0, earliest_failure_time)), 
                              FUN = createSurvGQMatrices, 
                              patient_age = patient_age,
                              baseline_hazard_knots=object$control$knots, 
                              baseline_hazard_ordSpline=object$control$ordSpline)
  
  #Now we have a function which creates XBeta_Wgamma etc. matrices
  Xbeta_Wgamma_h0BsGamma_Creator = function(betas_dre, betas_psa, gammas, Bs_gammas){
    obsDre_Xbeta = obsDre_X %*% betas_dre
    obsPsa_Xbeta = obsPsa_X %*% betas_psa
    obsPsa_cat_Xbeta = obsPsa_cat_X %*% betas_psa
    obs_exp_Wgamma = exp(obs_W %*% gammas)
    
    GQ_hazard_fixed_part=lapply(GQ_hazard_X_Z_W_h0, function(GQ_hazard_data){
      if(!is.null(GQ_hazard_data)){
        GQ_weights_baselinehazard = c(GQ_hazard_data$GQ_weights * exp(GQ_hazard_data$GQ_h0_matrix %*% Bs_gammas))
        
        GQ_psa_value_Xbeta = GQ_hazard_data$GQ_psa_value_X %*% betas_psa
        GQ_psa_slope_Xbeta = GQ_hazard_data$GQ_psa_slope_X_Z %*% betas_psa[4:7]
        GQ_dre_Xbeta = GQ_hazard_data$GQ_dre_X %*% betas_dre
        
        return(list(GQ_weights_baselinehazard = GQ_weights_baselinehazard,
                    GQ_psa_value_Xbeta = GQ_psa_value_Xbeta,
                    GQ_psa_slope_Xbeta = GQ_psa_slope_Xbeta,
                    GQ_dre_Xbeta = GQ_dre_Xbeta,
                    GQ_psa_value_Z = GQ_hazard_data$GQ_psa_value_Z,
                    GQ_psa_slope_Z = GQ_hazard_data$GQ_psa_slope_X_Z,
                    GQ_dre_Z = GQ_hazard_data$GQ_dre_Z))
      }else{
        return(NULL)
      }
    })
    
    return(list(obsDre_Xbeta=obsDre_Xbeta, obsPsa_Xbeta=obsPsa_Xbeta, obsPsa_cat_Xbeta=obsPsa_cat_Xbeta,
                obs_exp_Wgamma = obs_exp_Wgamma, GQ_hazard_fixed_part=GQ_hazard_fixed_part))
  }
  
  #############################################################################
  #1. First we find a proposal density for random effects using Empirical Bayes
  #############################################################################
  psa_cat_quad = lapply(lower_upper_psa_limits, FUN = function(limits_at_time){
    lapply(limits_at_time, FUN = function(limits){
      if(diff(limits)==Inf){
        return(limits)
      }else{
        return(getGaussianQuadWeightsPoints(limits))
      }
    })
  })
  
  log_numerator_bayesrule_data = list(obs_log2psaplus1 = patient_data_psa$log2psaplus1,
                                      obs_log2psaplus1_category = patient_data_psa_cat$log2psaplus1,
                                      obs_psa_cat_visitTimeYears = patient_data_psa_cat$visitTimeYears,
                                      obs_high_dre = patient_data_dre$high_dre,
                                      obsPsa_Z = obsPsa_Z, obsDre_Z = obsDre_Z,
                                      obsPsa_cat_Z = obsPsa_cat_Z,
                                      psa_cat_quad=psa_cat_quad,
                                      max_psa_cat = length(psa_cat_quad[[1]]),
                                      psaDist = psaDist,
                                      df = TdistDf,
                                      latest_survival_time=latest_survival_time, 
                                      earliest_failure_time=earliest_failure_time)
  
  log_numerator_bayesrule_data$inv_D = object$statistics$postMeans$inv_D
  log_numerator_bayesrule_data$sigma = object$statistics$postMeans$sigma2
  log_numerator_bayesrule_data$invSigmaMatrix = diag(nrow(patient_data_psa)) / object$statistics$postMeans$sigma2^2
  log_numerator_bayesrule_data$alphas = object$statistics$postMeans$alphas
  log_numerator_bayesrule_data$Xbeta_Wgamma_h0BsGamma = Xbeta_Wgamma_h0BsGamma_Creator(object$statistics$postMeans$betas1,
                                                                                       object$statistics$postMeans$betas2,
                                                                                       object$statistics$postMeans$gammas,
                                                                                       object$statistics$postMeans$Bs_gammas)
  
  optim_function <- function (b, data_list){-log_numerator_bayesrule(b, data_list)}
  gradient_function <- function (b, data_list){cd(b, optim_function, 
                                                  data = log_numerator_bayesrule_data, eps = 1e-03)}
  #par is start values
  for(optim_method in optim_methods){  
    empiricalbayes_b = try(optim(par = rep(0, ncol(object$statistics$postMeans$inv_D)),
                             fn = optim_function, 
                             gr = gradient_function, data = log_numerator_bayesrule_data, 
                             method = optim_method, hessian = TRUE,
                             control = list(maxit = 200,
                                            parscale = rep(0.1, ncol(object$statistics$postMeans$inv_D)))), silent = T)
    
    if(!inherits(empiricalbayes_b, 'try-error') && all(!is.nan(empiricalbayes_b$hessian))){
      break
    }
  }
  if(M==0 | inherits(empiricalbayes_b, 'try-error')){
    return(empiricalbayes_b)
  }
  
  #This has M elements
  postMCMC_theta_indices = sample(1:length(object$mcmc$sigma2), size = M)
  #Sample b from a proposal distribution
  #invVars_b <- opt$hessian / scale
  proposed_b = rmvt(n = M, mu=empiricalbayes_b$par, Sigma = scale * solve(empiricalbayes_b$hessian), df=4)
  log_pdf_proposed_b = dmvt(x=proposed_b, mu=empiricalbayes_b$par, Sigma = scale * solve(empiricalbayes_b$hessian), df=4, log=T)
  
  #Indices range between 1 and M+1
  current_b = empiricalbayes_b$par
  log_pdf_current_b = dmvt(x=current_b, mu=empiricalbayes_b$par, Sigma = scale * solve(empiricalbayes_b$hessian), df=4, log=T)
  
  accepted_b = matrix(NA, nrow = M, ncol = ncol(proposed_b))
  for(m in 1:M){
    theta_index = postMCMC_theta_indices[m]
    
    #First create the data for the new thetas
    log_numerator_bayesrule_data$inv_D = object$mcmc$inv_D[theta_index,,]
    log_numerator_bayesrule_data$sigma = object$mcmc$sigma2[theta_index]
    log_numerator_bayesrule_data$invSigmaMatrix = diag(nrow(patient_data_psa)) / object$mcmc$sigma2[theta_index]^2
    log_numerator_bayesrule_data$alphas = object$mcmc$alphas[theta_index,]
    log_numerator_bayesrule_data$Xbeta_Wgamma_h0BsGamma = Xbeta_Wgamma_h0BsGamma_Creator(object$mcmc$betas1[theta_index,],
                                                                                         object$mcmc$betas2[theta_index,],
                                                                                         object$mcmc$gammas[theta_index,],
                                                                                         object$mcmc$Bs_gammas[theta_index,])
    #log likelihood current value
    logl_data_current_b = log_numerator_bayesrule(b=current_b, data = log_numerator_bayesrule_data)
    
    #log likelihood next value
    logl_data_proposed_b = log_numerator_bayesrule(b=proposed_b[m,], data = log_numerator_bayesrule_data)
    
    acceptance_probability = min(exp(logl_data_proposed_b - logl_data_current_b + log_pdf_current_b - log_pdf_proposed_b[m]),1)
    if(!is.nan(acceptance_probability) & rbinom(n = 1, size = 1, prob = acceptance_probability)==1){
      accepted_b[m,] = proposed_b[m,]
      current_b = proposed_b[m,]
      log_pdf_current_b = log_pdf_proposed_b[m]
    }else{
      accepted_b[m,] = current_b
      #No changes in current_b, and log_pdf_current_b
    }
  }
  
  return(list(posterior_b=accepted_b, posterior_theta_mcmc_indices=postMCMC_theta_indices))
}

#the parameter addRandomError is for adding error term to the mean effect to 
#get a prediction interval
getExpectedFuture_Cat = function(object, patient_data, lower_upper_psa_limits=list(),
                                 latest_survival_time=0, earliest_failure_time=Inf, 
                                 survival_predict_times=NULL, Y_predict_times=NULL,
                                 psaDist = "normal", 
                                 TdistDf=3, M=200, rounding = 5){
  repeat{
    post_b_beta = try(get_b_fullBayes(object, patient_data, lower_upper_psa_limits,
                                latest_survival_time, earliest_failure_time,
                                scale = 1.6, psaDist, TdistDf=TdistDf, M), silent = T)
    if(!inherits(post_b_beta, 'try-error')){
      break
    }
    
    earliest_failure_time = earliest_failure_time + abs(jitter(0, amount = 0.05))
  }
  
  #M==0 check empirical bayes
  if(M==0){
    posterior_b = matrix(post_b_beta$par)
    mcmc_betas_dre = matrix(object$statistics$postMeans$betas1)
    mcmc_betas_psa = matrix(object$statistics$postMeans$betas2)
    mcmc_sigma_psa = object$statistics$postMeans$sigma2
    
    if(!is.null(survival_predict_times)){
      mcmc_alphas = matrix(object$statistics$postMeans$alphas)
      mcmc_Bs_gammas = matrix(object$statistics$postMeans$Bs_gammas)
      mcmc_gammas = matrix(object$statistics$postMeans$gammas)
    }
    M=1
  }else{
    posterior_b = t(post_b_beta$posterior_b)
    mcmc_betas_dre = t(object$mcmc$betas1[post_b_beta$posterior_theta_mcmc_indices,])
    mcmc_betas_psa = t(object$mcmc$betas2[post_b_beta$posterior_theta_mcmc_indices,])
    mcmc_sigma_psa = object$mcmc$sigma2[post_b_beta$posterior_theta_mcmc_indices,]
    
    if(!is.null(survival_predict_times)){
      mcmc_alphas = t(object$mcmc$alphas[post_b_beta$posterior_theta_mcmc_indices,])
      mcmc_Bs_gammas = t(object$mcmc$Bs_gammas[post_b_beta$posterior_theta_mcmc_indices,])
      mcmc_gammas = t(object$mcmc$gammas[post_b_beta$posterior_theta_mcmc_indices,])
    }
  }
  
  predicted_Y_prob = predicted_surv_prob = NULL
  
  if(!is.null(Y_predict_times) & length(Y_predict_times)>0){
    predicted_psa_means = psaXbetaZb(patient_data$Age[1], Y_predict_times, mcmc_betas_psa, posterior_b[3:7,])
    
    predicted_highdre_prob = plogis(dreLogOddsXbetaZb(patient_data$Age[1], Y_predict_times, mcmc_betas_dre, posterior_b[1:2,]))
    predicted_lowdre_prob = 1 - predicted_highdre_prob
    
    #Calculate PSA pnorm
    predicted_psacat_pnorm = lapply(1:NR_DISCRETIZED_PSA, function(psa_cat){
      ret = t(sapply(1:length(Y_predict_times), FUN = function(Ypt_index){
        psa_predict_time = Y_predict_times[Ypt_index]
        limits_at_time = lower_upper_psa_limits[[as.character(psa_predict_time)]]
        pnorm(q = limits_at_time[[psa_cat]][2], mean = predicted_psa_means[Ypt_index,], sd = mcmc_sigma_psa, lower.tail = T)
      }))
      if(M==1){
        ret = t(ret)
      }
      return(ret)
    })
    
    predicted_psacat_prob = vector("list", length = NR_DISCRETIZED_PSA)
    predicted_psacat_prob[[1]] = predicted_psacat_pnorm[[1]]
    predicted_psacat_prob[2:NR_DISCRETIZED_PSA] = lapply(2:NR_DISCRETIZED_PSA, function(x){
      predicted_psacat_pnorm[[x]] - predicted_psacat_pnorm[[x-1]]
    })
    
    predicted_Y_prob = vector("list", NR_DISCRETIZED_PSA * 2)
    predicted_Y_prob[1:NR_DISCRETIZED_PSA] = lapply(predicted_psacat_prob, function(psa_cat_prob){
      psa_cat_prob * predicted_lowdre_prob
    })
    predicted_Y_prob[(NR_DISCRETIZED_PSA+1):(2*NR_DISCRETIZED_PSA)] = lapply(predicted_psacat_prob, function(psa_cat_prob){
      psa_cat_prob * predicted_highdre_prob
    })
    
    if(length(Y_predict_times)>1){
      predicted_Y_prob = t(sapply(predicted_Y_prob, FUN = rowMeans))
    }else{
      predicted_Y_prob = matrix(sapply(predicted_Y_prob, FUN = rowMeans),
                                nrow=length(OUTCOME_CAT_NAMES), ncol = 1)
    }
    predicted_Y_prob = round(predicted_Y_prob, digits=rounding)
    rownames(predicted_Y_prob) = OUTCOME_CAT_NAMES
    colnames(predicted_Y_prob) = Y_predict_times
  }
  
  if(!is.null(survival_predict_times) & length(survival_predict_times)>0){
    survival_predict_times = survival_predict_times[survival_predict_times > latest_survival_time & survival_predict_times < earliest_failure_time]
    
    integration_time_pairs = rep(c(latest_survival_time, survival_predict_times, earliest_failure_time), each=2)
    integration_time_pairs = integration_time_pairs[!integration_time_pairs %in% c(0, Inf)]
    integration_time_pairs = c(0, integration_time_pairs[-length(integration_time_pairs)])
    
    integration_time_pairs = split(integration_time_pairs, rep(1:(length(integration_time_pairs)/2), each=2))
    
    patient_age = patient_data$Age[1]
    baseline_hazard_knots=object$control$knots
    baseline_hazard_ordSpline=object$control$ordSpline
    
    #here c(0,1) is just an arbitrary choice to check how many Quad points are used
    n_quad_points = length(getGaussianQuadWeightsPoints(c(0,1))$points)
    alpha_dre_matrix = matrix(mcmc_alphas[1,], nrow = n_quad_points, ncol=M, byrow = T)
    alpha_psa_matrix = matrix(mcmc_alphas[2,], nrow = n_quad_points, ncol=M, byrow = T)
    alpha_psa_slope_matrix = matrix(mcmc_alphas[3,], nrow = n_quad_points, ncol=M, byrow = T)
    W_gammas = model.matrix(survivalFormula, data = data.frame(Age = rep(patient_age,n_quad_points))) %*% mcmc_gammas
    
    cum_hazard_M = rep(0, M)
    surv_prob_M = list()
    for(i in 1:length(integration_time_pairs)){
      wp = getGaussianQuadWeightsPoints(integration_time_pairs[[i]])
      GQ_h0_matrix =  splineDesign(knots = baseline_hazard_knots, ord = baseline_hazard_ordSpline,
                                   x = wp$points, outer.ok = T)
      
      GQ_h0 = exp(GQ_h0_matrix %*% mcmc_Bs_gammas)
      fitted_psa_alpha = psaXbetaZb(patient_age, wp$points, mcmc_betas_psa, posterior_b[3:7,]) * alpha_psa_matrix
      fitted_psa_slope_alpha = psaSlopeXbetaZb(patient_age, wp$points, mcmc_betas_psa[4:7,], posterior_b[4:7,]) * alpha_psa_slope_matrix
      fitted_dre_alpha = dreLogOddsXbetaZb(patient_age, wp$points, mcmc_betas_dre, posterior_b[1:2,]) * alpha_dre_matrix
      total_hazard = GQ_h0 * exp(W_gammas + fitted_psa_alpha + fitted_psa_slope_alpha + fitted_dre_alpha)
      
      cum_hazard_M = cum_hazard_M + apply(total_hazard, 2, FUN = function(x){sum(x * wp$weights)})
      surv_prob_M[[i]] = exp(-cum_hazard_M)
    }
    
    if(latest_survival_time==0){
      denominator_left = rep(1, M)
    }else{
      denominator_left = surv_prob_M[[1]]
      surv_prob_M = surv_prob_M[2:length(surv_prob_M)]
    }
    
    if(earliest_failure_time==Inf){
      denominator_right = rep(0,M)
    }else{
      denominator_right = surv_prob_M[[length(surv_prob_M)]]
      surv_prob_M = surv_prob_M[-length(surv_prob_M)]
    }
    
    denominator = denominator_left - denominator_right
    if(M==1){
      predicted_surv_prob = matrix(sapply(surv_prob_M, function(x){(x-denominator_right)/denominator}), byrow = F)
    }else{
      predicted_surv_prob = t(sapply(surv_prob_M, function(x){(x-denominator_right)/denominator}))
    }
    rownames(predicted_surv_prob) = survival_predict_times
    
    predicted_surv_prob = apply(predicted_surv_prob,1, mean)
  }
  
  return(list(predicted_surv_prob=predicted_surv_prob,
              predicted_Y_prob=predicted_Y_prob))
}

environment(get_b_fullBayes) = asNamespace("JMbayes")
environment(getExpectedFuture_Cat) = asNamespace("JMbayes")
