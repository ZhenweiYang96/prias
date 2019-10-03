args = commandArgs(trailingOnly = T)

library(JMbayes)
library(survival)
library(splines)

cohort = args[1]
iteration = as.numeric(args[2])
seed = 2019 + iteration

load("Rdata/gap3/PRIAS_2019/motherdata.Rdata")
load(paste0("Rdata/gap3/PRIAS_2019/validation/recalibrated_prias_model/mvJoint_psa_recalib_prias",cohort, ".Rdata"))
load(paste0("Rdata/gap3/PRIAS_2019/validation/fitted_true_models/mvJoint_psa_",cohort, ".Rdata"))
source("src/clinical_gap3/prediction_only_psa.R")
source("src/clinical_gap3/validation/auc_brier/gof.R")

M = 750
max_thoriz = round(reclassification_df$time_10pat_risk_set[reclassification_df$Cohort==cohort])

print(paste("Cohort: ", cohort, "-- seed: ", seed, "-- max_thoriz: ", max_thoriz))

t_horizs = seq(1, max_thoriz, 0.5)
set.seed(seed)

longdata = get(paste0("longdata_", cohort))
if(iteration > 1){
  print("bootstrapped iteration")  
  
  bs_oldpids = sample(unique(longdata$P_ID), replace = T)
  
  longdata = do.call('rbind', lapply(1:length(bs_oldpids), FUN = function(j){
    tempds = longdata[longdata$P_ID == bs_oldpids[j],]
    tempds$P_ID = j
    return(tempds)
  }))
}else{
  print('original dataset iteration')
}

pe_list = vector("list", length(t_horizs))

for(k in 1:length(pe_list)){
  print(paste("Calculations started for thoriz:", t_horizs[k]))
  
  print('Calculating PE')
  pe_list[[k]] = goodness_of_fit(pred_model = prias_model_recalib,
                                   orig_model = mvJoint_psa_time_scaled,
                                   newdata = longdata,
                                   T_start = t_horizs[k]-1,
                                   T_horiz = t_horizs[k],
                                   M = M)
  
  #save(pe_list, file = paste0("Rdata/gap3/PRIAS_2019/pe/",cohort,"_",seed,".Rdata"))
}
