library(ggplot2)
library(JMbayes)
library(splines)

load("Rdata/gap3/PRIAS_2019/validation/fitted_true_models/mvJoint_psa_MUSIC.Rdata")
cohort_model = mvJoint_psa_time_scaled
load("Rdata/gap3/PRIAS_2019/mvJoint_psa_time_scaled_light.Rdata")
prias_model = mvJoint_psa_time_scaled

source("src/clinical_gap3/prediction_only_psa.R")

longdata = cohort_model$model_info$mvglmer_components$data
longdata.id = longdata[!duplicated(longdata$P_ID),]
longdata.id$right_cens_time = longdata.id$latest_survival_time
longdata.id$right_cens_time[longdata.id$earliest_failure_time!=Inf] = (0.5 * (longdata.id$latest_survival_time + longdata.id$earliest_failure_time))[longdata.id$earliest_failure_time!=Inf]
longdata.id$reclassification = longdata.id$earliest_failure_time!=Inf

nr_knots = 6
cohort_model_ordSpline = cohort_model$control$ordSpline

tt = longdata.id$right_cens_time[longdata.id$reclassification==T]
pp = quantile(tt, c(0.05, 0.95), names = FALSE)
kn = tail(head(seq(pp[1L], pp[2L], length.out = nr_knots), -1L), -1L)
kn <- kn[kn < max(longdata.id$right_cens_time)]
st = outer(tt/2, sk + 1)
cohort_model_knots = sort(c(rep(range(longdata.id$right_cens_time, st), cohort_model_ordSpline), kn))

nr_bs_gammas = dim(splineDesign(cohort_model_knots, seq(0.1,0.9, length.out = 10),cohort_model_ordSpline))[2]

#Now we optimize the log likelihood of time to event data
#this time we estimate the baseline hazard parameters

#Now we create for each patient all data that is needed before optim
optim_data_list = lapply(split(longdata, f = longdata$P_ID), function(pat_data){
  set.seed(2019)
  
  latest_survival_time = pat_data$latest_survival_time[1]
  earliest_failure_time = pat_data$earliest_failure_time[1]
  patient_age = pat_data$age[1]
  
  M=500
  #Random effects for this patient using PRIAS model
  post_b_beta = get_b_fullBayes(object = prias_model,
                                patient_data = pat_data,
                                latest_survival_time = latest_survival_time,
                                earliest_failure_time = earliest_failure_time,
                                psaDist = "Tdist", M = M)
  
  posterior_b = t(post_b_beta$posterior_b)
  mcmc_betas_psa = t(prias_model$mcmc$betas1[post_b_beta$posterior_theta_mcmc_indices,])
  mcmc_alphas = t(prias_model$mcmc$alphas[post_b_beta$posterior_theta_mcmc_indices,])
  mcmc_gammas = t(prias_model$mcmc$gammas[post_b_beta$posterior_theta_mcmc_indices,])
  
  createSurvGQMatrices = function(lower_upper_limit){
    if(sum(lower_upper_limit) != Inf & diff(lower_upper_limit)!=0){
      wp = getGaussianQuadWeightsPoints(lower_upper_limit)
      survival_predict_times = wp$points
      
      bh_matrix =  splineDesign(knots = cohort_model_knots, ord = cohort_model_ordSpline,
                                x = survival_predict_times, outer.ok = T)
      
      W_gammas = model.matrix(survivalFormula, data = data.frame(age = rep(patient_age, length(survival_predict_times)))) %*% mcmc_gammas
      
      alpha_psa_matrix = matrix(mcmc_alphas[1,], nrow = length(survival_predict_times), ncol=M, byrow = T)
      alpha_psa_slope_matrix = matrix(mcmc_alphas[2,], nrow = length(survival_predict_times), ncol=M, byrow = T)
      
      fitted_psa_alpha = psaXbetaZb(patient_age, survival_predict_times, mcmc_betas_psa, posterior_b) * alpha_psa_matrix
      fitted_psa_slope_alpha = psaSlopeXbetaZb(patient_age, survival_predict_times, mcmc_betas_psa[-c(1:2),], posterior_b[-1,]) * alpha_psa_slope_matrix
      linear_predictor = W_gammas + fitted_psa_alpha + fitted_psa_slope_alpha
      
      weighted_mean_exp_lp = wp$weights * exp(rowMeans(linear_predictor))
      
      return(list(bh_matrix = bh_matrix,
                  weighted_mean_exp_lp = weighted_mean_exp_lp))
    }else{
      return(NULL)
    }
  }
  
  GQ_hazard_partial = lapply(list("latest_survival_time_data"=c(0, latest_survival_time),
                                  "earliest_failure_time_data"=c(0, earliest_failure_time)), 
                             FUN = createSurvGQMatrices)
  
  return(list(latest_survival_time=latest_survival_time,
              earliest_failure_time=earliest_failure_time,
              GQ_hazard_partial=GQ_hazard_partial))
})

cumulative_hazard = function(bs_gammas, GQ_hazard_partial){
  baseline_hazard = exp(GQ_hazard_partial$bh_matrix %*% bs_gammas)
  cumulativeHazard = sum(baseline_hazard * GQ_hazard_partial$weighted_mean_exp_lp)
  return(cumulativeHazard)
}

optim_function = function(bs_gammas, optim_data){
  
  log_likelikehoods = sapply(optim_data, FUN = function(data_list){
    cum_hazard_latest_survival_time = if(data_list$latest_survival_time==0){
      0
    }else{
      cumulative_hazard(bs_gammas, data_list$GQ_hazard_partial$latest_survival_time_data)
    }
    
    cum_hazard_earliest_failure_time = if(data_list$earliest_failure_time==Inf){
      Inf
    }else{
      cumulative_hazard(bs_gammas, data_list$GQ_hazard_partial$earliest_failure_time_data)
    }
    
    return(log(exp(-cum_hazard_latest_survival_time) - exp(-cum_hazard_earliest_failure_time)))
  })
  
  return(-sum(log_likelikehoods, na.rm = T))
}

gradient_function <- function (bs_gammas, optim_data){
  JMbayes:::cd(bs_gammas, optim_function, 
               optim_data = optim_data, eps = 1e-03)
}
for(optim_method in c("L-BFGS-B", "BFGS", "CG")){  
  mle_bs_gammas = try(optim(par = rep(0, nr_bs_gammas),
                            fn = optim_function, 
                            gr = gradient_function, optim_data = optim_data_list, 
                            method = optim_method, hessian = TRUE,
                            control = list(maxit = 500,
                                           parscale = rep(0.1, nr_bs_gammas))), silent = T)
  
  if(!inherits(mle_bs_gammas, 'try-error') && all(!is.nan(mle_bs_gammas$hessian))){
    break
  }
}

prias_model_recalib = prias_model
prias_model_recalib$statistics$postMeans$Bs_gammas = mle_bs_gammas$par
prias_model_recalib$control$knots = cohort_model_knots
prias_model_recalib$control$ordSpline = cohort_model_ordSpline
prias_model_recalib$mcmc$Bs_gammas = MASS::mvrnorm(n=nrow(prias_model_recalib$mcmc$Bs_gammas),
                                                   mu = mle_bs_gammas$par,
                                                   Sigma = 1.6*solve(mle_bs_gammas$hessian))

save(prias_model_recalib, file="Rdata/gap3/PRIAS_2019/validation/recalibrated_prias_model/mvJoint_psa_recalib_priasMUSIC.Rdata")