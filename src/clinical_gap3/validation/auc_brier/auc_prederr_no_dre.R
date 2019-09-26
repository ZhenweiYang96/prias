library(JMbayes)
library(survival)
library(splines)

load("Rdata/gap3/PRIAS_2019/cleandata.Rdata")
load("Rdata/gap3/PRIAS_2019/mvJoint_psa_time_scaled.Rdata")
source("src/clinical_gap3/auc_brier/auc_mod_prias.R")
source("src/clinical_gap3/auc_brier/prederr_mod_prias.R")

nr_bootstrap = 26

M = 750
t_horizs = seq(1, 10, 0.5)
center_name = "PRIAS"
set.seed(2019)

#This is done to help in bootstrap dataset generation
dataset_list = vector("list", length = nr_bootstrap)
#First item is original item
dataset_list[[1]] = prias_long_final

dataset_list[2:nr_bootstrap] = lapply(2:nr_bootstrap, FUN = function(x){
  bs_oldpids = sample(unique(prias_long_final$P_ID), replace = T)
  
  newds = do.call('rbind', lapply(1:length(bs_oldpids), FUN = function(j){
    tempds = prias_long_final[prias_long_final$P_ID == bs_oldpids[j],]
    tempds$P_ID = j
    return(tempds)
  }))
  
  return(newds)
})

auc_prederr_bs = vector("list", nr_bootstrap)

for(j in 1:nr_bootstrap){
  print(paste("Calculations started for bootstrapped sample:", j))
  auc_prederr_bs[[j]]$auc_list = vector("list", length(t_horizs))
  auc_prederr_bs[[j]]$prederr_list = vector("list", length(t_horizs))
  
  for(k in 1:length(t_horizs)){
    auc_prederr_bs[[j]]$auc_list[[k]] = aucJM.mvJMbayes_mod(mvJoint_psa_time_scaled,
                                                            newdata = dataset_list[[j]],
                                                            Tstart = t_horizs[k]-1,
                                                            Thoriz = t_horizs[k],
                                                            idVar = "P_ID", M = M)
    
    print(paste0('AUC for ', t_horizs[k]-1,"--", t_horizs[k], ": ", auc_prederr_bs[[j]]$auc_list[[k]]$auc))
    save(auc_prederr_bs, file = paste0("Rdata/gap3/PRIAS_2019/auc_prederr/",center_name,".Rdata"))
    
    auc_prederr_bs[[j]]$prederr_list[[k]] = prederrJM.mvJMbayes_mod(mvJoint_psa_time_scaled, 
                                                                    newdata = dataset_list[[j]], 
                                                                    Tstart = t_horizs[k]-1, 
                                                                    Thoriz = t_horizs[k], 
                                                                    idVar = "P_ID", M = M)
    
    print(paste0('Prederr for ', t_horizs[k]-1,"--", t_horizs[k], ": ", auc_prederr_bs[[j]]$prederr_list[[k]]$prederr))
    save(auc_prederr_bs, file = paste0("Rdata/gap3/PRIAS_2019/auc_prederr/", center_name,".Rdata"))
  }
}
