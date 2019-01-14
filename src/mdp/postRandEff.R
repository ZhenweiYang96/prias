psaXbetaZb = function(Age, visitTimeYears, betas_psa, b_psa){
  psaXbeta = function(Age, visitTimeYears, betas_psa){
    fixed_psaFormula = ~ 1 +I(Age - 70) +  I((Age - 70)^2) + ns(visitTimeYears, knots=c(0.1, 0.7, 4), Boundary.knots=c(0, 5.42))
    df = data.frame(Age, visitTimeYears)
    c(model.matrix(fixed_psaFormula, df) %*% betas_psa)
  }
  
  random_psa_formula = ~ 1 + ns(visitTimeYears, knots=c(0.1, 0.7, 4), Boundary.knots=c(0, 5.42))
  random_psa_df = model.matrix(random_psa_formula, data.frame(visitTimeYears))
  
  psaXbeta(Age, visitTimeYears, betas_psa) + random_psa_df %*% b_psa
}

psaSlopeXbetaZb = function(Age, visitTimeYears, betas_psa_slope, b_psa_slope){
  fixed_random_psaSlopeFormula = ~ 0 + dns(visitTimeYears, knots=c(0.1, 0.7, 4), Boundary.knots=c(0, 5.42))
  
  X_Z_matrix = model.matrix(fixed_random_psaSlopeFormula, data.frame(Age, visitTimeYears))
  
  X_Z_matrix %*% betas_psa_slope + X_Z_matrix %*% b_psa_slope
}

dreLogOddsXbetaZb = function(Age, visitTimeYears, betas_dre, b_dre){
  dreLogOddsXbeta = function(Age, visitTimeYears, betas_dre){
    fixed_dreFormula = ~ 1 + I(Age - 70) +  I((Age - 70)^2) + visitTimeYears
    df = data.frame(Age, visitTimeYears)
    
    c(model.matrix(fixed_dreFormula, df) %*% betas_dre)
  }
  
  random_dre_formula = ~ 1 + visitTimeYears
  random_dre_df = model.matrix(random_dre_formula, data.frame(visitTimeYears))
  
  dreLogOddsXbeta(Age, visitTimeYears, betas_dre) + random_dre_df %*% b_dre
}

surv_prob = function(b, obs_exp_Wgamma, GQ_hazard_data, alphas){
  GQ_hazard_true_dreLogOdds = GQ_hazard_data$GQ_dre_Xbeta + GQ_hazard_data$GQ_dre_Z %*% b[1:2]  
  GQ_hazard_true_psa = GQ_hazard_data$GQ_psa_value_Xbeta + GQ_hazard_data$GQ_psa_value_Z %*% b[3:7]
  GQ_hazard_true_psaSlope = GQ_hazard_data$GQ_psa_slope_Xbeta + GQ_hazard_data$GQ_psa_slope_Z %*% b[4:7]
  
  GQ_hazard_longitudinal_part = exp(GQ_hazard_true_dreLogOdds * alphas[1] + GQ_hazard_true_psa * alphas[2] + GQ_hazard_true_psaSlope * alphas[3])
  
  cumulativeHazard = sum(GQ_hazard_data$GQ_weights_baselinehazard * obs_exp_Wgamma * GQ_hazard_longitudinal_part)
  return(exp(-cumulativeHazard))
}

log_p_T = function(b, timeLowerLimit, timeUpperLimit, obs_exp_Wgamma, GQ_hazard_X_Z_W_h0_theta, alphas){
  surv_prob_timeLowerLimit = if(timeLowerLimit==0){
    1
  }else{
    surv_prob(b, obs_exp_Wgamma, GQ_hazard_X_Z_W_h0_theta$lowerLimitData, alphas)
  }
  
  surv_prob_timeUpperLimit = if(timeUpperLimit==Inf){
    0
  }else{
    surv_prob(b, obs_exp_Wgamma, GQ_hazard_X_Z_W_h0_theta$upperLimitData, alphas)
  }
  
  return(log(surv_prob_timeLowerLimit - surv_prob_timeUpperLimit))
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
  
  #Time to event contribution
  time_to_event_contrib = log_p_T(b, data$timeLowerLimit, data$timeUpperLimit, data$Xbeta_Wgamma_h0BsGamma$obs_exp_Wgamma,
                                  data$Xbeta_Wgamma_h0BsGamma$GQ_hazard_X_Z_W_h0_theta, data$alphas)
  
  return(randomEff_prior_Contrib + dre_contrib + psa_contrib + time_to_event_contrib)
}

#When M=0 it becomes empirical bayes
get_b_fullBayes = function(object, patient_data, timeLowerLimit, timeUpperLimit=Inf,
                           scale = 1.6, psaDist = "normal", TdistDf=3, M=200, seed = 2019){
  
  set.seed(seed)
  
  fixed_psaFormula = ~ 1 +I(Age - 70) +  I((Age - 70)^2) + ns(visitTimeYears, knots=c(0.1, 0.7, 4), Boundary.knots=c(0, 5.42))
  random_psaFormula = ~ 1 + ns(visitTimeYears, knots=c(0.1, 0.7, 4), Boundary.knots=c(0, 5.42))
  fixed_random_psaSlopeFormula = ~ 0 + dns(visitTimeYears, knots=c(0.1, 0.7, 4), Boundary.knots=c(0, 5.42))
  
  fixed_dreFormula = ~ 1 + I(Age - 70) +  I((Age - 70)^2) + visitTimeYears
  random_dreFormula = ~ 1 + visitTimeYears
  
  survivalFormula = ~ 0 + I(Age - 70) + I((Age - 70)^2)
  
  #Obs data X and Z matrices for longitudinal part
  patient_data_dre = patient_data[!is.na(patient_data$high_dre), c("Age",  "visitTimeYears", "high_dre")]
  patient_data_psa = patient_data[!is.na(patient_data$log2psaplus1), c("Age",  "visitTimeYears", "log2psaplus1")]
  patient_age = patient_data$Age[1]
  
  obsPsa_X = model.matrix(fixed_psaFormula, patient_data_psa)
  obsPsa_Z = model.matrix(random_psaFormula, patient_data_psa)
  obsDre_X = model.matrix(fixed_dreFormula, patient_data_dre)
  obsDre_Z = model.matrix(random_dreFormula, patient_data_dre)
  
  #obs W, and baseline hazard matrix for time to event part
  #For time to event contribution: basically calculate hazard and multiply by weights
  GQsurv_nodes_weights = if(object$control$GQsurv == "GaussKronrod"){
    gaussKronrod()
  }else{ 
    gaussLegendre(control$GQsurv.k)
  }
  
  obs_W = model.matrix(survivalFormula, data = data.frame(Age = patient_age))
  GQ_hazard_X_Z_W_h0 = lapply(c("lowerLimitData"=timeLowerLimit, "upperLimitData"=timeUpperLimit), 
                              FUN = function(time){
                                if(!time %in% c(0,Inf)){
                                  GQ_scaled_weights = GQsurv_nodes_weights$wk * time/2
                                  GQ_scaled_visitTimeYears = c(outer(time/2, GQsurv_nodes_weights$sk + 1))
                                  
                                  GQ_h0_matrix =  splineDesign(knots = object$control$knots, ord = object$control$ordSpline,
                                                               x = GQ_scaled_visitTimeYears, outer.ok = T)
                                  
                                  GQ_patient_df = data.frame(Age = patient_age, visitTimeYears = GQ_scaled_visitTimeYears)
                                  
                                  GQ_psa_value_X = model.matrix(fixed_psaFormula, GQ_patient_df)
                                  GQ_psa_value_Z = model.matrix(random_psaFormula, GQ_patient_df)
                                  GQ_psa_slope_X_Z = model.matrix(fixed_random_psaSlopeFormula, GQ_patient_df)
                                  GQ_dre_X = model.matrix(fixed_dreFormula, GQ_patient_df)
                                  GQ_dre_Z = model.matrix(random_dreFormula, GQ_patient_df)
                                  
                                  return(list(GQ_scaled_weights = GQ_scaled_weights,
                                              GQ_h0_matrix = GQ_h0_matrix,
                                              GQ_psa_value_X = GQ_psa_value_X,
                                              GQ_psa_value_Z = GQ_psa_value_Z,
                                              GQ_psa_slope_X_Z = GQ_psa_slope_X_Z,
                                              GQ_dre_X = GQ_dre_X,
                                              GQ_dre_Z = GQ_dre_Z))
                                }else{
                                  return(NULL)
                                }
                              })
  
  #Now we have a function which creates XBeta_Wgamma etc. matrices
  Xbeta_Wgamma_h0BsGamma_Creator = function(betas_dre, betas_psa, gammas, Bs_gammas){
    obsDre_Xbeta = obsDre_X %*% betas_dre
    obsPsa_Xbeta = obsPsa_X %*% betas_psa
    obs_exp_Wgamma = exp(obs_W %*% gammas)
    
    GQ_hazard_X_Z_W_h0_theta=lapply(GQ_hazard_X_Z_W_h0, function(GQ_hazard_data){
      if(!is.null(GQ_hazard_data)){
        GQ_weights_baselinehazard = c(GQ_hazard_data$GQ_scaled_weights * exp(GQ_hazard_data$GQ_h0_matrix %*% Bs_gammas))
        
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
    
    return(list(obsDre_Xbeta=obsDre_Xbeta, obsPsa_Xbeta=obsPsa_Xbeta,
                obs_exp_Wgamma = obs_exp_Wgamma, GQ_hazard_X_Z_W_h0_theta=GQ_hazard_X_Z_W_h0_theta))
  }
  
  #############################################################################
  #1. First we find a proposal density for random effects using Empirical Bayes
  #############################################################################
  log_numerator_bayesrule_data = list(obs_log2psaplus1 = patient_data_psa$log2psaplus1,
                                      obs_high_dre = patient_data_dre$high_dre,
                                      obsPsa_Z=obsPsa_Z, obsDre_Z=obsDre_Z,
                                      psaDist = psaDist,
                                      df = TdistDf,
                                      timeLowerLimit=timeLowerLimit, 
                                      timeUpperLimit=timeUpperLimit)
  
  log_numerator_bayesrule_data$inv_D = object$statistics$postMeans$inv_D
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
  empiricalbayes_b = optim(par = rep(0, ncol(object$statistics$postMeans$inv_D)),
                           fn = optim_function, 
                           gr = gradient_function, data = log_numerator_bayesrule_data, 
                           method = "BFGS", hessian = TRUE, 
                           control = list(maxit = 200, 
                                          parscale = rep(0.1, ncol(object$statistics$postMeans$inv_D))))
  
  if(M==0){
    return(empiricalbayes_b)
  }
  
  #Sample b from a proposal distribution
  #invVars_b <- opt$hessian / scale
  proposed_b = rmvt(n = M, mu=empiricalbayes_b$par, Sigma = scale * solve(empiricalbayes_b$hessian), df=4)
  log_pdf_proposed_b = dmvt(x=proposed_b, mu=empiricalbayes_b$par, Sigma = scale * solve(empiricalbayes_b$hessian), df=4, log=T)
  
  #This has M elements
  postMCMC_theta_indices = sample(1:length(object$mcmc$sigma2), size = M)
  
  #Indices range between 1 and M+1
  current_b = empiricalbayes_b$par
  log_pdf_current_b = dmvt(x=current_b, mu=empiricalbayes_b$par, Sigma = scale * solve(empiricalbayes_b$hessian), df=4, log=T)
  
  accepted_b = matrix(NA, nrow = M, ncol = ncol(proposed_b))
  for(m in 1:M){
    theta_index = postMCMC_theta_indices[m]
    
    #First create the data for the new thetas
    log_numerator_bayesrule_data$inv_D = object$mcmc$inv_D[theta_index,,]
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

predictLongitudinalOutcome = function(object, patient_data, timeLowerLimit, timeUpperLimit=Inf, visitTimeYears,
                                      psaDist = "normal", TdistDf=3, M=200, seed = 2019){
  
  post_b_beta = get_b_fullBayes(object, patient_data, timeLowerLimit, timeUpperLimit,
                                scale = 1.6, psaDist, TdistDf=TdistDf, M, seed)  
  
  posterior_b = post_b_beta$posterior_b
  posterior_theta_mcmc_indices = post_b_beta$posterior_theta_mcmc_indices
  
  mcmc_betas_dre = object$mcmc$betas1
  mcmc_betas_psa = object$mcmc$betas2
  
  predicted_psa = sapply(1:M, FUN = function(m){
    psaXbetaZb(patient_data$Age[1], visitTimeYears, 
               mcmc_betas_psa[posterior_theta_mcmc_indices[m],], posterior_b[m,3:7])
  })
  
  predicted_psa_slope = sapply(1:M, FUN = function(m){
    psaSlopeXbetaZb(patient_data$Age[1], visitTimeYears, 
               mcmc_betas_psa[posterior_theta_mcmc_indices[m],4:7], posterior_b[m,4:7])
  })
  
  predicted_dreLogOdds = sapply(1:M, FUN = function(m){
    dreLogOddsXbetaZb(patient_data$Age[1], visitTimeYears, 
               mcmc_betas_dre[posterior_theta_mcmc_indices[m],], posterior_b[m,1:2])
  })
  
  return(list(predicted_psa=predicted_psa, predicted_psa_slope=predicted_psa_slope, predicted_dreLogOdds=predicted_dreLogOdds))
}
