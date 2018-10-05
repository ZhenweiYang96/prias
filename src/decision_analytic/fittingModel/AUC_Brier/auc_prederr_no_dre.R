load("Rdata/decision_analytic/cleandata.Rdata")
source("src/decision_analytic/load_lib.R")
source("src/decision_analytic/fittingModel/AUC_Brier/auc_mod_prias.R")
source("src/decision_analytic/fittingModel/AUC_Brier/prederr_mod_prias.R")

maxCores = 6

totalSizeBootstrap = 20

list.prias.id = vector("list", totalSizeBootstrap)
list.prias_long = vector("list", totalSizeBootstrap)

auc_no_dre = vector("list", totalSizeBootstrap)
prederr_no_dre = vector("list", totalSizeBootstrap)

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
    
    print(paste("Dataset generated", i))
    
    new.survModel_no_dre = survreg(Surv(progression_time_start, progression_time_end, type = "interval2") ~ 
                                  I(Age - 70) +  I((Age - 70)^2), data = list.prias.id[[i]], model = TRUE)
    
    print(paste("survival model fitted", i))
    
    psa_data_set_i =  list.prias_long[[i]][!is.na(list.prias_long[[i]]$psa),]
    
    new.mvglm_no_dre = try(mvglmer(list(log2psaplus1 ~ I(Age - 70) +  I((Age - 70)^2) +
                                       ns(visitTimeYears, knots=c(0.1, 0.7, 4), Boundary.knots=c(0, 5.42)) + 
                                       (ns(visitTimeYears, knots=c(0.1, 0.7, 4), Boundary.knots=c(0, 5.42))|P_ID)),
                                data=psa_data_set_i, families = list(gaussian)))
    
    if(inherits(new.mvglm_no_dre, "try-error")==FALSE){
      break
    }else{
      print("Trying again after mvglmer error")
    }
  }
  
  print(paste("mvglmer fitted", i))
  print(paste("Seed used:", lastSeed))
  
  seeds = c(seeds, lastSeed)
  
  new.mvJoint_no_dre = mvJointModelBayes(new.mvglm_no_dre, new.survModel_no_dre, 
                                      timeVar = "visitTimeYears", Formulas = list("log2psaplus1" = "value",
                                                                                  "log2psaplus1" = list(fixed = ~ 0 + dns(visitTimeYears, knots=c(0.1, 0.7, 4), Boundary.knots=c(0, 5.42)),
                                                                                                        random=~0 + dns(visitTimeYears, knots=c(0.1, 0.7, 4), Boundary.knots=c(0, 5.42)),
                                                                                                        indFixed = 4:7, indRandom=2:5, name = "slope")),
                                      control=list(n_cores = maxCores))
  
  print(paste("JM fitted PSA, PSA velocity", i))
  
  times = 1:9 + 0.0001
  deltaT = 1.0101
  
  auc_no_dre[[i]] = vector("list", length(times))
  prederr_no_dre[[i]] = vector("list", length(times))
  
  j=1
  for(Tstart in times){
    auc_no_dre[[i]][[j]] = aucJM.mvJMbayes_mod(new.mvJoint_no_dre, list.prias_long[[i]],
                             Tstart = Tstart, Thoriz=Tstart + deltaT, 
                             idVar="P_ID", n_cores = maxCores)
    
    prederr_no_dre[[i]][[j]] = prederrJM.mvJMbayes_mod(new.mvJoint_no_dre, list.prias_long[[i]],
                                                       Tstart = Tstart, Thoriz=Tstart + deltaT, 
                                                       idVar="P_ID")
    
    print(paste("AUC no_dre, time=", Tstart, i))
    
    j = j + 1
    
    auc_prederr_no_dre = list("auc_no_dre" = auc_no_dre, "prederr_no_dre"=prederr_no_dre)
    save(auc_prederr_no_dre, file="Rdata/decision_analytic/AUC_Brier/auc_prederr_no_dre.Rdata")
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