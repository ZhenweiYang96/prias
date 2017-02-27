library(survival)

#FIrst: Choose only those patients who have Gleason <=6 at baseline
condition_survival_gleason1 = prias_long$visit_number==1 & prias_long$gleason %in% c(7,8,9,10)
#no need to using unique as only 1 visit is used in the condition
prias_long_survival_eligible = prias_long[!(prias_long$P_ID %in% droplevels(prias_long[condition_survival_gleason1, "P_ID"])) ,]

prias_long_survival_eligible$P_ID = droplevels(prias_long_survival_eligible$P_ID)
length(unique(prias_long_survival_eligible$P_ID))

#Second: Choose only those patients who have atleast one gleason score
condition_survival_gleason2pid = droplevels(unique(prias_long_survival_eligible[!is.na(prias_long_survival_eligible$gleason),]$P_ID))
prias_long_survival_eligible = prias_long_survival_eligible[prias_long_survival_eligible$P_ID %in% condition_survival_gleason2pid,]
prias_long_survival_eligible$P_ID = droplevels(prias_long_survival_eligible$P_ID)
length(unique(prias_long_survival_eligible$P_ID))

#Third:  drop all observations after the first time Gleason score is > 7
prias.id = prias_long_survival_eligible[!duplicated(prias_long_survival_eligible$P_ID),]
pid_progressed = droplevels(unique(prias_long_survival_eligible[prias_long_survival_eligible$gleason %in% c(7,8,9,10),]$P_ID))

prias.id$progressed = ifelse(prias.id$P_ID %in% pid_progressed, 1, 0)
prias.id$progression_time = by(prias_long_survival_eligible, prias_long_survival_eligible$P_ID, function(prias_long_i){
  if(any(prias_long_i$gleason %in% c(7,8,9,10))){
    prias_long_i$visitTimeYears[which(prias_long_i$gleason %in% c(7,8,9,10))[1]]  
  }else{
    max(prias_long_i$visitTimeYears, na.rm = T)
  }
})

prias_long_survival_eligible$progression_time = unlist(lapply(prias.id$P_ID, function(pid){
  rowCount = nrow(prias_long_survival_eligible[prias_long_survival_eligible$P_ID == pid,])
  rep(prias.id[prias.id$P_ID==pid,"progression_time"][1], rowCount)
}))

prias_long_survival_eligible = prias_long_survival_eligible[prias_long_survival_eligible$visitTimeYears<=prias_long_survival_eligible$progression_time,]

##############################
prias.id$AgeGroup = sapply(prias.id$Age, function(age){
  if(age<=67){
    return("<=67")
  }else if(age>73){
    return(">73")
  }else{
    return(">67 & <=73")
  }
})

kmfit = survfit(Surv(progression_time, progressed)~AgeGroup, conf.type="log-log", data=prias.id)
survminer::ggsurvplot(kmfit, risk.table = T,break.time.by = 1, xlab = "Time(years)", ylim = c(0.5,1), conf.int = T)

coxModel = coxph(Surv(progression_time, progressed) ~ I(Age - 70) +  I((Age - 70)^2), 
                 data=prias.id, x = T, model = T)
