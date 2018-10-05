load("Rdata/decision_analytic/cleandata.Rdata")
source("src/decision_analytic/load_lib.R")
source("src/decision_analytic/fittingModel/AUC_Brier/auc_mod_prias.R")
source("src/decision_analytic/fittingModel/AUC_Brier/prederr_mod_prias.R")

maxCores = 6

#Step 1: AUC and Prederr on Real prias data with real prias fitted model object
auc_0_1 = aucJM.mvJMbayes_mod(mvJoint_dre_psa_dre_value, newdata = prias_long, 
                    Tstart = 0.0001, Thoriz = 1.0001, idVar = "P_ID")
auc_1_2 = aucJM.mvJMbayes_mod(mvJoint_dre_psa_dre_value, newdata = prias_long, 
                    Tstart = 1.0001, Thoriz = 2.0001, idVar = "P_ID")
auc_2_3 = aucJM.mvJMbayes_mod(mvJoint_dre_psa_dre_value, newdata = prias_long, 
                    Tstart = 2.0001, Thoriz = 3.0001, idVar = "P_ID")
auc_3_4 = aucJM.mvJMbayes_mod(mvJoint_dre_psa_dre_value, newdata = prias_long, 
                    Tstart = 3.0001, Thoriz = 4.0001, idVar = "P_ID")
auc_4_5 = aucJM.mvJMbayes_mod(mvJoint_dre_psa_dre_value, newdata = prias_long, 
                    Tstart = 4.0001, Thoriz = 5.0001, idVar = "P_ID")

save.image(file="Rdata/decision_analytic/AUC_Brier/auc_prederr_t_1.Rdata")

prederr_0_1 = prederrJM.mvJMbayes_mod(mvJoint_dre_psa_dre_value, newdata = prias_long, 
                              Tstart = 0.0001, Thoriz = 1.0001, idVar = "P_ID")
prederr_1_2 = prederrJM.mvJMbayes_mod(mvJoint_dre_psa_dre_value, newdata = prias_long, 
                              Tstart = 1.0001, Thoriz = 2.0001, idVar = "P_ID")
prederr_2_3 = prederrJM.mvJMbayes_mod(mvJoint_dre_psa_dre_value, newdata = prias_long, 
                              Tstart = 2.0001, Thoriz = 3.0001, idVar = "P_ID")
prederr_3_4 = prederrJM.mvJMbayes_mod(mvJoint_dre_psa_dre_value, newdata = prias_long, 
                              Tstart = 3.0001, Thoriz = 4.0001, idVar = "P_ID")
prederr_4_5 = prederrJM.mvJMbayes_mod(mvJoint_dre_psa_dre_value, newdata = prias_long, 
                              Tstart = 4.0001, Thoriz = 5.0001, idVar = "P_ID")

save.image(file="Rdata/decision_analytic/AUC_Brier/auc_prederr_t_1.Rdata")



totalSizeBootstrap = 20

list.prias.id = vector("list", totalSizeBootstrap)
list.prias_long = vector("list", totalSizeBootstrap)

aucAll = vector("list", totalSizeBootstrap)
prederrAll = vector("list", totalSizeBootstrap)

getNextSeed = function(lastSeed){
  lastSeed + 1
}

seeds = c()
lastSeed = 2018
for(i in 1:totalSizeBootstrap){
  while(T){
    lastSeed = getNextSeed(lastSeed)
    set.seed(lastSeed)
    
    sampledRows = sort(sample(1:nrow(prias.id), replace = T))
    old_PIDs = prias.id$P_ID[sampledRows]
    
    list.prias.id[[i]] = prias.id[sampledRows,]
    list.prias.id[[i]]$P_ID = 1:length(sampledRows)
    
    list.prias_long[[i]] = do.call(rbind, args = lapply(list.prias.id[[i]]$P_ID, function(newId){
      temp = prias_long[prias_long$P_ID == old_PIDs[newId],]
      temp$P_ID = newId
      temp$progressed = list.prias.id[[i]]$progressed[list.prias.id[[i]]$P_ID == newId]
      temp$progression_time_start = list.prias.id[[i]]$progression_time_start[list.prias.id[[i]]$P_ID == newId]
      temp$progression_time_end = list.prias.id[[i]]$progression_time_end[list.prias.id[[i]]$P_ID == newId]
      
      return(temp)
    }))
    
    list.prias_long[[i]]$high_dre = ifelse(list.prias_long[[i]]$dre=="T1c", 0, 1)
    
    print(paste("Dataset generated", i))
    
    new.survModel_all = survreg(Surv(progression_time_start, progression_time_end, type = "interval2") ~ 
                              I(Age - 70) +  I((Age - 70)^2), data = list.prias.id[[i]], model = TRUE)
      
    print(paste("survival model fitted", i))
    
    dre_psa_data_set_i =  list.prias_long[[i]][!(is.na(list.prias_long[[i]]$dre) & is.na(list.prias_long[[i]]$psa)),]
    
    new.mvglm_all = try(mvglmer(list(high_dre~I(Age - 70) +  I((Age - 70)^2) + visitTimeYears  + 
                   (visitTimeYears|P_ID),
                 
                 log2psaplus1 ~ I(Age - 70) +  I((Age - 70)^2) +
                   ns(visitTimeYears, knots=c(0.1, 0.7, 4), Boundary.knots=c(0, 5.42)) + 
                   (ns(visitTimeYears, knots=c(0.1, 0.7, 4), Boundary.knots=c(0, 5.42))|P_ID)),
            
            data=dre_psa_data_set_i, families = list(binomial, gaussian), engine = "STAN",
            control = list(n.iter=1000)))
    
    if(inherits(new.mvglm_all, "try-error")==FALSE){
      break
    }else{
      print("Trying again after mvglmer error")
    }
  }
  
  print(paste("mvglmer fitted", i))
  print(paste("Seed used:", lastSeed))
  
  seeds = c(seeds, lastSeed)
  
  forms_dre_value = list("high_dre" = "value",
                         "log2psaplus1" = "value",
                         "log2psaplus1" = list(fixed = ~ 0 + dns(visitTimeYears, knots=c(0.1, 0.7, 4), Boundary.knots=c(0, 5.42)),
                                               random=~0 + dns(visitTimeYears, knots=c(0.1, 0.7, 4), Boundary.knots=c(0, 5.42)),
                                               indFixed = 4:7, indRandom=2:5, name = "slope"))
  
  new.mvJoint_all = mvJointModelBayes(new.mvglm_all, new.survModel_all, 
                                  timeVar = "visitTimeYears", Formulas = forms_dre_value,
                                  control=list(n_cores = maxCores))
  
  print(paste("JM fitted DRE, PSA, PSA velocity", i))
  
  times = 1:9 + 0.0001
  deltaT = 1.0101
  
  aucAll[[i]] = vector("list", length(times))
  prederrAll[[i]] = vector("list", length(times))
  
  j=1
  for(Tstart in times){
    aucAll[[i]][[j]] = aucJM.mvJMbayes_mod(new.mvJoint_all, list.prias_long[[i]],
                                    Tstart = Tstart, Thoriz=Tstart + deltaT, 
                                 idVar="P_ID", n_cores = maxCores)
    
    prederrAll[[i]][[j]] = prederrJM.mvJMbayes_mod(new.mvJoint_all, list.prias_long[[i]],
                                                       Tstart = Tstart, Thoriz=Tstart + deltaT, 
                                                       idVar="P_ID")
    
    
    print(paste("AUC All, time=", Tstart, i))
    
    j = j + 1
    
    auc_prederr_all = list("auc_all" = aucAll, "prederr_all"=prederrAll)
    save(auc_prederr_all, file="Rdata/decision_analytic/AUC_Brier/auc_prederr_all.Rdata")
  }
}

# aucTdBoth = t(sapply(aucBootStrapped$both, function(x){
#   c(x[[1]]$auc, x[[2]]$auc, x[[3]]$auc)
# }))
# 
# aucTdVal = t(sapply(aucBootStrapped$val, function(x){
#   c(x[[1]]$auc, x[[2]]$auc, x[[3]]$auc)
# }))
# 
# aucTdBothMean = apply(aucTdBoth, MARGIN = 2, mean)
# aucTdBothLow = apply(aucTdBoth, MARGIN = 2, quantile, probs=0.025)
# aucTdBothUp = apply(aucTdBoth, MARGIN = 2, quantile, probs=0.975)
# 
# aucTdValMean = apply(aucTdVal, MARGIN = 2, mean)
# aucTdValLow = apply(aucTdVal, MARGIN = 2, quantile, probs=0.025)
# aucTdValUp = apply(aucTdVal, MARGIN = 2, quantile, probs=0.975)
# 
# plotDf = data.frame(time=rep(c(1,2,3),2), aucMean = c(aucTdBothMean, aucTdValMean),
#                     aucLow = c(aucTdBothLow, aucTdValLow), aucUp = c(aucTdBothUp, aucTdValUp),
#                     model=rep(c("with both log (PSA) value and velocity", 
#                                 "with only log (PSA) value"), each=3))
# 
# ggplot(data=plotDf) + geom_line(aes(x=time, y=aucMean, linetype=model)) +
#   geom_ribbon(aes(x=time, ymin=aucLow, ymax=aucUp, color=model), fill = "grey", alpha=0.4) + 
#   scale_linetype_manual(values=c("twodash", "solid")) +
#   theme(legend.title=element_blank(), legend.position = "top") + 
#   xlab("Time (years)") + ylab("Area under the curve (AUC)")
# 