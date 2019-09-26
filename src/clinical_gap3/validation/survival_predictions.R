#We compare 3 types of predictions
#linear predictor and calibration slope
#survival probability: weak and strong 

library(JMbayes)
library(splines)
library(survival)
library(interval) 

load("Rdata/gap3/PRIAS_2019/validation/fitted_true_models/mvJoint_psa_MUSIC.Rdata")
cohort_model = mvJoint_psa_time_scaled
load("Rdata/gap3/PRIAS_2019/mvJoint_psa_time_scaled_light.Rdata")
prias_model = mvJoint_psa_time_scaled
rm(mvJoint_psa_time_scaled)
load("Rdata/gap3/PRIAS_2019/validation/recalibrated_prias_model/mvJoint_psa_recalib_priasMUSIC.Rdata")

source("src/clinical_gap3/prediction_only_psa.R")

longdata = cohort_model$model_info$mvglmer_components$data
longdata.id = longdata[!duplicated(longdata$P_ID),]
longdata.id$right_cens_time = longdata.id$latest_survival_time
longdata.id$right_cens_time[longdata.id$earliest_failure_time!=Inf] = (0.5 * (longdata.id$latest_survival_time + longdata.id$earliest_failure_time))[longdata.id$earliest_failure_time!=Inf]
longdata.id$reclassification = longdata.id$earliest_failure_time!=Inf

calib_pred_times = seq(0, 10, 0.1)

models = list(cohort_model, prias_model, prias_model_recalib)
cumrisk_models = lapply(models, FUN = function(model){
  cumrisk = lapply(split(longdata, f = longdata$P_ID), function(pat_data){
    set.seed(2019)
    
    latest_survival_time = pat_data$latest_survival_time[1]
    earliest_failure_time = pat_data$earliest_failure_time[1]
    patient_age = pat_data$age[1]
    
    M=500
    
    surv_probs = rep(NA, length(calib_pred_times))
    surv_probs[calib_pred_times<=latest_survival_time] = 1
    surv_probs[calib_pred_times>=earliest_failure_time] = 0
    
    if(nrow(pat_data) > 1 & any(is.na(surv_probs))){
      temp = getExpectedFutureOutcomes(object = model,
                                       patient_data = pat_data,
                                       latest_survival_time = latest_survival_time,
                                       earliest_failure_time = earliest_failure_time,
                                       survival_predict_times = calib_pred_times,
                                       psaDist = "Tdist", M = M)
      if(any(is.na(surv_probs))){
        surv_probs[calib_pred_times > latest_survival_time &
                     calib_pred_times < earliest_failure_time] = rowMeans(temp$predicted_surv_prob)
      }
    }
    return(1-surv_probs)
  })
  cumrisk = do.call('cbind', cumrisk)
  rownames(cumrisk) = calib_pred_times
  
  return(cumrisk)
})
names(cumrisk_models) = c("cohort_model", "prias_model", "prias_model_recalib")

save(cumrisk_models, file="Rdata/gap3/PRIAS_2019/validation/predicted_risk_comparisons/MUSIC.Rdata")
