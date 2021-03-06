args = commandArgs(trailingOnly = T)

library(JMbayes)
library(survival)
library(splines)

cohort = args[1]
iteration = as.numeric(args[2])
seed = 2019 + iteration

load("Rdata/gap3/PRIAS_2019/motherdata.Rdata")
load("Rdata/gap3/PRIAS_2019/mvJoint_psa_time_scaled.Rdata")
source("src/clinical_gap3/validation/auc_brier/auc_mod_prias.R")

#For recalibrated
load(paste0("Rdata/gap3/PRIAS_2019/validation/recalibrated_prias_model/mvJoint_psa_recalib_prias",
            cohort,".Rdata"))
mvJoint_psa_time_scaled$mcmc$Bs_gammas = prias_model_recalib$mcmc$Bs_gammas
mvJoint_psa_time_scaled$statistics$postMeans$Bs_gammas = prias_model_recalib$statistics$postMeans$Bs_gammas
mvJoint_psa_time_scaled$control$knots = prias_model_recalib$control$knots
#

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

auc_list = vector("list", length(t_horizs))

for(k in 1:length(auc_list)){
  print(paste("Calculations started for thoriz:", t_horizs[k]))
  
  auc_list[[k]] = aucJM.mvJMbayes_mod(mvJoint_psa_time_scaled,
                                      newdata = longdata,
                                      Tstart = t_horizs[k]-1,
                                      Thoriz = t_horizs[k],
                                      idVar = "P_ID", M = M)
  
  print(paste0('AUC for ', t_horizs[k]-1,"--", t_horizs[k], ": ", auc_list[[k]]$auc))
  save(auc_list, file = paste0("Rdata/gap3/PRIAS_2019/auc/",cohort,"_",seed,".Rdata"))
}
