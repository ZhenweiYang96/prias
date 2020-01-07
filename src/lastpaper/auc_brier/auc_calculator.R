args = commandArgs(trailingOnly = T)

library(JMbayes)
library(survival)
library(splines)

iteration = as.numeric(args[1])
seed = 2019 + iteration

source("src/lastpaper/auc_brier/auc_mod_prias.R")
load("Rdata/gap3/PRIAS_2019/cleandata.Rdata")
load("Rdata/lastpaper/fitted_model/mvJoint_dre_psa_2knots_quad_age.Rdata")

M = 400
max_thoriz = 6

print(paste("seed: ", seed, "-- max_thoriz: ", max_thoriz))

t_horizs = seq(1, max_thoriz, 0.5)
set.seed(seed)

if(iteration > 1){
  print("bootstrapped iteration")  
  
  bs_oldpids = sample(unique(prias_long_final$P_ID), replace = T)
  
  prias_long_final = do.call('rbind', lapply(1:length(bs_oldpids), FUN = function(j){
    tempds = prias_long_final[prias_long_final$P_ID == bs_oldpids[j],]
    tempds$P_ID = j
    return(tempds)
  }))
}else{
  print('original dataset iteration')
}

auc_list = vector("list", length(t_horizs))

for(k in 1:length(auc_list)){
  print(paste("Calculations started for thoriz:", t_horizs[k]))
  
  auc_list[[k]] = aucJM.mvJMbayes_mod(mvJoint_dre_psa_2knots_quad_age,
                                      newdata = prias_long_final,
                                      Tstart = t_horizs[k]-1,
                                      Thoriz = t_horizs[k],
                                      idVar = "P_ID", M = M)
  
  print(paste0('AUC for ', t_horizs[k]-1,"--", t_horizs[k], ": ", auc_list[[k]]$auc))
  save(auc_list, file = paste0("Rdata/lastpaper/auc_pe/auc_",seed,".Rdata"))
}
