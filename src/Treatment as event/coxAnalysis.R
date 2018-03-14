library(survival)

long_ds = gleason_data_set
long_ds = dre_data_set
long_ds = psa_data_set

long_ds$P_ID = droplevels(long_ds$P_ID)

pid_cumsum = cumsum(table(long_ds$P_ID))
first_row_index_eachsub = c(0, pid_cumsum[-length(pid_cumsum)]) + 1
prias.id = long_ds[first_row_index_eachsub, c(1:6,9,12)]
prias.id$event_indicator = ifelse(prias.id$event_type %in% c("Treatment-Progression", "Died-Progression"),yes = 1, no = 0)
sum(prias.id$event_indicator)
write.csv(prias.id, file = "prias_id.csv", row.names = F)

prias.id$AgeGroup = sapply(prias.id$Age, function(age){
  if(age<=67){
    return("<=67")
  }else if(age>73){
    return(">73")
  }else{
    return(">67 & <=73")
  }
})

kmfit = survfit(Surv(year_discontinued, event_indicator)~AgeGroup, conf.type="log-log", data=prias.id)
survMisc:::autoplot.survfit(kmfit) + xlab("Time (years)")
survminer::ggsurvplot(kmfit, risk.table = T,break.time.by = 1, xlab = "Time(years)", ylim = c(0.5,1), conf.int = T)


#Ignore patient 957, no gleason for him

#This one gives a negative variance when calcuating from the inverse of the hessian
coxModel = coxph(Surv(year_discontinued, event_indicator) ~ Age + I(Age^2), 
                 data=prias.id, x = T, model = T)

coxModel = coxph(Surv(year_discontinued, event_indicator) ~ I(Age - 70) +  I((Age - 70)^2), 
                 data=prias.id, x = T, model = T)

coxModel = coxph(Surv(year_discontinued, event_indicator) ~ poly(Age,2)[,1] + poly(Age,2)[,2], 
                 data=prias.id, x = T, model = T)
anova(coxModel)

#Cox snell residuals
coxsnellres=prias.id$event_indicator-resid(coxModel,type="martingale")
fitres=survfit(coxph(Surv(coxsnellres,prias.id$event_indicator)~1,method='breslow'),type='aalen')
plot(fitres$time,-log(fitres$surv),type='s',xlab='Cox-Snell Residuals', 
     ylab='Estimated Cumulative Hazard Function')
abline(0,1,col='red',lty=2)

