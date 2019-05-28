library(JMbayes)
library(survival)
library(splines)

load("Rdata/prias/mvJoint_dre_psa_dre_value.Rdata")
source("src/auc_brier/auc_mod_prias.R")
source("src/auc_brier/prederr_mod_prias.R")
load('Rdata/longdata.Rdata')

large_centers = names(sort(by(longdata$center,
                              data = longdata$P_ID, 
                              function(x){length(unique(x))}), 
                           decreasing = T))[2:10]
large_centers = c(large_centers, "prias_orig", "PRIAS")

set.seed(2019)
M = 500
t_starts = 0:9
for(center_name in large_centers){

  load(paste0("Rdata/data/longdata_", center_name,".Rdata"))
  longds = get(paste0("longdata_", center_name))
  
  print(paste0("Doing for center: ", center_name))
  
  longds$progression_time_start = longds$latest_survival_time
  longds$progression_time_end = longds$earliest_failure_time
  
  auc_list = vector("list", length(t_starts))
  prederr_list = vector("list", length(t_starts))
  for(t_start in t_starts){
    auc_list[[t_start + 1]] = aucJM.mvJMbayes_mod(mvJoint_dre_psa_dre_value,
                                                  newdata = longds,
                                                   Tstart = t_start,
                                                   Thoriz = t_start + 1,
                                                   idVar = "P_ID", M = M)

    print(paste0('AUC for ', t_start,"--", t_start + 1, ": ", auc_list[[t_start + 1]]$auc))
    save(auc_list, file = paste0("Rdata/auc/auc_list_",center_name,".Rdata"))
    
    prederr_list[[t_start + 1]] = prederrJM.mvJMbayes_mod(mvJoint_dre_psa_dre_value, 
                                                  newdata = longds, 
                                                  Tstart = t_start, 
                                                  Thoriz = t_start + 1, 
                                                  idVar = "P_ID", M = M)
    
    print(paste0('Prederr for ', t_start,"--", t_start+1, ": ", prederr_list[[t_start + 1]]$prederr))
    save(prederr_list, file = paste0("Rdata/prederr/prederr_list_", center_name,".Rdata"))
  }
}
