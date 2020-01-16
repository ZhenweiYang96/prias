library(icenReg)

load("Rdata/gap3/PRIAS_2019/motherdata.Rdata")
source("src/clinical_gap3/validation/auc_brier/gof_intercept_only_model.R")

for(cohort in cohortnames){
  
  longdata_orig = get(paste0("longdata_", cohort))
  model = ic_sp(Surv(latest_survival_time, earliest_failure_time, type = 'interval2') ~ 1,
                data = longdata_orig)
  
  for(iteration in 1:30){
    seed = 2019 + iteration
    
    # plot(model, axes=F)
    # axis(side = 1, at = seq(1,15,0.5))
    # axis(side = 2, at = seq(0,1,0.1))
    # box()
    
    max_thoriz = round(reclassification_df$time_10pat_risk_set[reclassification_df$Cohort==cohort])
    
    print(paste("Cohort: ", cohort, "-- seed: ", seed, "-- max_thoriz: ", max_thoriz))
    
    t_horizs = seq(1, max_thoriz, 0.5)
    set.seed(seed)
    
    longdata = longdata_orig
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
      pe_list[[k]] = goodness_of_fit_intercept_only(model,
                                                    newdata = longdata,
                                                    T_start = t_horizs[k]-1,
                                                    T_horiz = t_horizs[k])
      
      save(pe_list, file = paste0("Rdata/gap3/PRIAS_2019/validation/pe_intercept_only_model/",cohort,"_",seed,".Rdata"))
    }
  }
}