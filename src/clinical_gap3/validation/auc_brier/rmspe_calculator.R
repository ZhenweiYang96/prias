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
max_auc_time = round(reclassification_df$time_10pat_risk_set[reclassification_df$Cohort==cohort])
t_horizs = seq(1, max_auc_time, 0.5)
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
}

print(paste("Calculations started for bootstrapped sample:", j))
prederr_list = vector("list", length(t_horizs))

for(k in 1:length(auc_list)){
  auc_list[[k]] = aucJM.mvJMbayes_mod(prederrJM.mvJMbayes_mod,
                                      newdata = longdata,
                                      Tstart = t_horizs[k]-1,
                                      Thoriz = t_horizs[k],
                                      idVar = "P_ID", M = M)
  
  print(paste0('Prederr for ', t_horizs[k]-1,"--", t_horizs[k], ": ", prederr_list[[k]]$prederr))
  save(auc_list, file = paste0("Rdata/gap3/PRIAS_2019/auc/",cohort,"_",seed,".Rdata"))
}


for(k in 1:length(t_horizs)){
  auc_prederr_bs[[j]]$prederr_list[[k]] = prederrJM.mvJMbayes_mod(mvJoint_psa_time_scaled, 
                                                                  newdata = dataset_list[[j]], 
                                                                  Tstart = t_horizs[k]-1, 
                                                                  Thoriz = t_horizs[k], 
                                                                  idVar = "P_ID", M = M)
  
  
  save(auc_prederr_bs, file = paste0("Rdata/gap3/PRIAS_2019/auc_prederr/", center_name,".Rdata"))
}
