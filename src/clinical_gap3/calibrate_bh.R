cohortname = "Toronto"

library(JMbayes)
library(splines)
library(survival)
library(interval)

#load("~/Desktop/ErasmusMC_datasets/PRIAS-2019/rocresults/rocresults.Rdata")
load("~/Data/ErasmusMC_datasets/PRIAS-2019/rocresults.Rdata")
load("Rdata/gap3/PRIAS_2019/npmle_all.Rdata")

longdata.id = get(paste0("longdata_", cohortname, ".id"))
longdata = get(paste0("longdata_", cohortname))

longdata.id$right_cens_time = longdata.id$latest_survival_time
longdata.id$right_cens_time[longdata.id$earliest_failure_time!=Inf] = (0.5 * (longdata.id$latest_survival_time + longdata.id$earliest_failure_time))[longdata.id$earliest_failure_time!=Inf]
longdata.id$reclassification = longdata.id$earliest_failure_time!=Inf
longdata$age = longdata$Age
longdata$year_visit = longdata$visitTimeYears
longdata$right_cens_time = longdata$latest_survival_time
longdata$right_cens_time[longdata$earliest_failure_time!=Inf] = (0.5 * (longdata$latest_survival_time + longdata$earliest_failure_time))[longdata$earliest_failure_time!=Inf]
longdata$reclassification = longdata$earliest_failure_time!=Inf

npmle=npmle_all[[cohortname]]
npmle_time_points = as.numeric(cbind(c(0,0), npmle$intmap))
npmle_cumrisk = c(0, as.numeric(rep(c(0, cumsum(npmle$pf)), each=2)))[1:length(npmle_time_points)]

rm(list = setdiff(ls(), c("longdata.id", "longdata", "cohortname", "npmle_time_points", "npmle_cumrisk")))

load("Rdata/gap3/PRIAS_2019/gap3cohorts/mvJoint_psa_Toronto.Rdata")
source("src/clinical_gap3/prediction_only_psa.R")
load(paste0("Rdata/gap3/PRIAS_2019/gap3cohorts/pred_pers_cum_risk_prias_in_gap3/personalized_cumrisk_jm_",cohortname, ".Rdata"))

#First get hazard per person using the Toronto model
cumrisk_jm = lapply(longdata.id$P_ID, function(id){
  print(id)
  latest_survival_time = longdata.id$latest_survival_time[longdata.id$P_ID==id]
  earliest_failure_time = longdata.id$earliest_failure_time[longdata.id$P_ID==id]

  surv_probs = rep(NA, length(altman_calib_pred_times))
  surv_probs[altman_calib_pred_times<=latest_survival_time] = 1
  surv_probs[altman_calib_pred_times>=earliest_failure_time] = 0

  if(!id %in% names(which(table(longdata$P_ID) == 1))){
    temp = getExpectedFutureOutcomes(object = mvJoint_psa_time_scaled,
                                     patient_data = longdata[longdata$P_ID==id,],
                                     latest_survival_time = latest_survival_time,
                                     earliest_failure_time = earliest_failure_time,
                                     survival_predict_times = altman_calib_pred_times,
                                     psaDist = "Tdist", M = 500)
    if(any(is.na(surv_probs))){
      surv_probs[altman_calib_pred_times > latest_survival_time &
                   altman_calib_pred_times < earliest_failure_time] = rowMeans(temp$predicted_surv_prob)
    }
  }
  return(1-surv_probs)
})
