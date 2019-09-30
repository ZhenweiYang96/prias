args = commandArgs(trailingOnly = T)

library(JMbayes)
library(survival)
library(splines)

cohort = args[1]
iteration = as.numeric(args[2])
seed = 2019 + iteration

load("Rdata/gap3/PRIAS_2019/motherdata.Rdata")
load("Rdata/gap3/PRIAS_2019/mvJoint_psa_time_scaled.Rdata")
source("src/clinical_gap3/validation/auc_brier/prederr_mod_prias.R")

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

mspe_list = vector("list", length(t_horizs))
mape_list = vector("list", length(t_horizs))

for(k in 1:length(mspe_list)){
  print(paste("Calculations started for thoriz:", t_horizs[k]))
  
  print('Calculating MSPE')
  mspe_list[[k]] = prederrJM.mvJMbayes_mod(mvJoint_psa_time_scaled,
                                      newdata = longdata,
                                      Tstart = t_horizs[k]-1,
                                      Thoriz = t_horizs[k],
                                      idVar = "P_ID", M = M)
  
  print(paste0('MSPE for ', t_horizs[k]-1,"--", t_horizs[k], ": ", mspe_list[[k]]$prederr))
  save(mspe_list, file = paste0("Rdata/gap3/PRIAS_2019/mspe/",cohort,"_",seed,".Rdata"))
  
  print('Calculating MAPE')
  mape_list[[k]] = prederrJM.mvJMbayes_mod(mvJoint_psa_time_scaled,
                                           newdata = longdata,
                                           Tstart = t_horizs[k]-1,
                                           Thoriz = t_horizs[k],
                                           idVar = "P_ID", M = M, lossFun='absolute')
  
  print(paste0('MAPE for ', t_horizs[k]-1,"--", t_horizs[k], ": ", mape_list[[k]]$prederr))
  save(mape_list, file = paste0("Rdata/gap3/PRIAS_2019/mape/",cohort,"_",seed,".Rdata"))
}
