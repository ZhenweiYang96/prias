load("Rdata/decision_analytic/DRE_PSA/mvJoint_dre_psa_dre_value.Rdata")
source("src/decision_analytic/load_lib.R")
source("src/decision_analytic/fittingModel/AUC_Brier/auc_mod_prias.R")
source("src/decision_analytic/fittingModel/AUC_Brier/prederr_mod_prias.R")

maxCores = 4

fileNames = list.files("Rdata/decision_analytic/AUC_Brier/all/", full.names = T)

for(fileName in fileNames){
  load(fileName)
    
  times = 0:4 + 0.0001
  deltaT = 1
  
  #bb means bootstrapped data, bootstrapped model. bo means bootstrapped data, original model
  aucAll_bb = vector("list", length(times))
  prederrAll_bb = vector("list", length(times))
  
  aucAll_bo = vector("list", length(times))
  prederrAll_bo = vector("list", length(times))
  
  j=1
  for(Tstart in times){
    aucAll_bb[[j]] = aucJM.mvJMbayes_mod(bootstrapdata$mvJoint, bootstrapdata$data_long,
                                           Tstart = Tstart, Thoriz=Tstart + deltaT, 
                                           idVar="P_ID", n_cores = maxCores)
    
    aucAll_bo[[j]] = aucJM.mvJMbayes_mod(mvJoint_dre_psa_dre_value, bootstrapdata$data_long,
                                         Tstart = Tstart, Thoriz=Tstart + deltaT, 
                                         idVar="P_ID", n_cores = maxCores)
    
    prederrAll_bb[[j]] = prederrJM.mvJMbayes_mod(bootstrapdata$mvJoint, bootstrapdata$data_long,
                                                   Tstart = Tstart, Thoriz=Tstart + deltaT, 
                                                   idVar="P_ID")
    
    prederrAll_bo[[j]] = prederrJM.mvJMbayes_mod(mvJoint_dre_psa_dre_value, bootstrapdata$data_long,
                                                 Tstart = Tstart, Thoriz=Tstart + deltaT, 
                                                 idVar="P_ID")
    
    bootstrapdata$aucAll_bb = aucAll_bb
    bootstrapdata$prederrAll_bb = prederrAll_bb
    
    bootstrapdata$aucAll_bo = aucAll_bo
    bootstrapdata$prederrAll_bo = prederrAll_bo
    
    save(bootstrapdata, file=fileName)
    j = j + 1
  }
}