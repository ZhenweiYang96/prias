args = commandArgs(trailingOnly = T)

library(JMbayes)
library(survival)
library(splines)

iteration = as.numeric(args[1])
seed = 2019 + iteration

load("Rdata/gap3/PRIAS_2019/cleandata.Rdata")
load("Rdata/lastpaper/fitted_model/mvJoint_dre_psa_2knots_quad_age.Rdata")
source("src/lastpaper/prediction_psa_dre.R")
source("src/lastpaper/auc_brier/gof.R")

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

pe_list = vector("list", length(t_horizs))

for(k in 1:length(pe_list)){
  print(paste("Calculations started for thoriz:", t_horizs[k]))
  
  print('Calculating PE')
  pe_list[[k]] = goodness_of_fit(pred_model = mvJoint_dre_psa_2knots_quad_age,
                                 orig_model = mvJoint_dre_psa_2knots_quad_age,
                                 newdata = prias_long_final,
                                 T_start = t_horizs[k]-1,
                                 T_horiz = t_horizs[k],
                                 M = M)
  
  save(pe_list, file = paste0("Rdata/lastpaper/auc_pe/pe_",seed,".Rdata"))
}
