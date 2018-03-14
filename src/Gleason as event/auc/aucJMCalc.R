list.prias.id = vector("list", 15)
list.prias_long = vector("list", 15)

load("Rdata/Gleason as event/auc/aucBootStrapped.Rdata")

#aucTdBoth = vector("list", 15)
#aucTdVal = vector("list", 15)
aucTdBoth = aucBootStrapped$both
aucTdVal = aucBootStrapped$val

seeds = c(1001, 1002)

getNextSeed = function(lastSeed){
  lastSeed + 1
}

lastSeed = 3310
for(i in 2:2){
  while(T){
    lastSeed = getNextSeed(lastSeed)
    set.seed(lastSeed)
    
    sampledRows = sort(sample(1:nrow(prias.id.rightCens), replace = T))
    old_PIDs = prias.id.rightCens$P_ID[sampledRows]
    
    list.prias.id[[i]] = prias.id.rightCens[sampledRows,]
    list.prias.id[[i]]$P_ID = 1:length(sampledRows)
    
    list.prias_long[[i]] = do.call(rbind, args = lapply(list.prias.id[[i]]$P_ID, function(newId){
      temp = prias_long[prias_long$P_ID == old_PIDs[newId],]
      temp$P_ID = newId
      temp$progressed = list.prias.id[[i]]$progressed[list.prias.id[[i]]$P_ID == newId]
      temp$progression_time = list.prias.id[[i]]$progression_time[list.prias.id[[i]]$P_ID == newId]
      
      return(temp)
    }))
    
    print(paste("Dataset generated", i))
    
    new.survModel_rightCens = coxph(Surv(progression_time, progressed) ~ I(Age - 70) +  I((Age - 70)^2), 
                                    data=list.prias.id[[i]], x = T, model = T)
    
    print(paste("Cox model fitted", i))
    
    new.mvglm_psa_spline_pt1pt54_pt1 = try(mvglmer(list(
      log2psa ~  I(Age - 70) +  I((Age - 70)^2) + 
        ns(visitTimeYears, knots=c(0.1, 0.5, 4), Boundary.knots=c(0, 7)) + 
        (ns(visitTimeYears, knots=c(0.1), Boundary.knots=c(0, 7))|P_ID)), 
      data=list.prias_long[[i]][!is.na(list.prias_long[[i]]$log2psa),], 
      families = list(gaussian)))
    
    if(inherits(new.mvglm_psa_spline_pt1pt54_pt1, "try-error")==FALSE){
      break
    }else{
      print("Trying again after mvglmer error")
    }
  }
  print(paste("mvglmer fitted", i))
  print(paste("Seed used:", lastSeed))
  
  seeds = c(seeds, lastSeed)
  
  new.mvJoint_psa_spline_pt1pt54_pt1_tdboth_rc = mvJointModelBayes(new.mvglm_psa_spline_pt1pt54_pt1, new.survModel_rightCens, timeVar = "visitTimeYears", 
                                                                   Formulas = list("log2psa" = "value",
                                                                                   "log2psa" = list(fixed = ~ 0 + dns(visitTimeYears, knots=c(0.1, 0.5, 4), Boundary.knots=c(0, 7)),
                                                                                                    random=~0 + dns(visitTimeYears, knots=c(0.1), Boundary.knots=c(0, 7)),
                                                                                                    indFixed = 4:7, indRandom=2:3, name = "slope")))
  
  print(paste("TD both JM fitted", i))
  
  new.mvJoint_psa_spline_pt1pt54_pt1_tdval_rc = mvJointModelBayes(new.mvglm_psa_spline_pt1pt54_pt1, new.survModel_rightCens, 
                                                                  timeVar = "visitTimeYears")
  
  print(paste("TD val JM fitted", i))
  
  aucTdBoth[[i]] = vector("list", 3)
  aucTdVal[[i]] = vector("list", 3)
  
  j=1
  for(Tstart in c(1.0001, 2.0001, 3.0001)){
    aucTdBoth[[i]][[j]] = aucJM_mod(new.mvJoint_psa_spline_pt1pt54_pt1_tdboth_rc,
                                    list.prias_long[[i]],
                                    Tstart = Tstart, Thoriz=Tstart + 1.0101, idVar="P_ID")
    
    print(paste("AUC JM TD both, time=", Tstart, i))
    
    aucTdVal[[i]][[j]] = aucJM_mod(new.mvJoint_psa_spline_pt1pt54_pt1_tdval_rc,
                                   list.prias_long[[i]],
                                   Tstart = Tstart, Thoriz=Tstart + 1.0101, idVar="P_ID")
    
    print(paste("AUC JM TD val, time=", Tstart, i))
    
    j = j + 1
  }
  
  aucBootStrapped = list("both" = aucTdBoth, "val"=aucTdVal)
  save(aucBootStrapped, file="Rdata/aucBootStrapped.Rdata")
}

aucTdBoth = t(sapply(aucBootStrapped$both, function(x){
  c(x[[1]]$auc, x[[2]]$auc, x[[3]]$auc)
}))

aucTdVal = t(sapply(aucBootStrapped$val, function(x){
  c(x[[1]]$auc, x[[2]]$auc, x[[3]]$auc)
}))

aucTdBothMean = apply(aucTdBoth, MARGIN = 2, mean)
aucTdBothLow = apply(aucTdBoth, MARGIN = 2, quantile, probs=0.025)
aucTdBothUp = apply(aucTdBoth, MARGIN = 2, quantile, probs=0.975)

aucTdValMean = apply(aucTdVal, MARGIN = 2, mean)
aucTdValLow = apply(aucTdVal, MARGIN = 2, quantile, probs=0.025)
aucTdValUp = apply(aucTdVal, MARGIN = 2, quantile, probs=0.975)

plotDf = data.frame(time=rep(c(1,2,3),2), aucMean = c(aucTdBothMean, aucTdValMean),
                    aucLow = c(aucTdBothLow, aucTdValLow), aucUp = c(aucTdBothUp, aucTdValUp),
                    model=rep(c("with both log (PSA) value and velocity", 
                                "with only log (PSA) value"), each=3))

ggplot(data=plotDf) + geom_line(aes(x=time, y=aucMean, linetype=model)) +
  geom_ribbon(aes(x=time, ymin=aucLow, ymax=aucUp, color=model), fill = "grey", alpha=0.4) + 
  scale_linetype_manual(values=c("twodash", "solid")) +
  theme(legend.title=element_blank(), legend.position = "top") + 
  xlab("Time (years)") + ylab("Area under the curve (AUC)")
