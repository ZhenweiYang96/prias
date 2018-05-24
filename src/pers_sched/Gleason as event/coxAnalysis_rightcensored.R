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
