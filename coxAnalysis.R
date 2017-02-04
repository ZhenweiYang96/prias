library(survival)

pid_cumsum = cumsum(table(prias_long$P_ID))
first_row_index_eachsub = c(0, pid_cumsum[-length(pid_cumsum)]) + 1
prias.id = prias_long[first_row_index_eachsub, c(1:6,9,12)]
prias.id$event_indicator = ifelse(prias.id$event_type %in% c("Treatment-Progression", "Died-Progression"),yes = 1, no = 0)
write.csv(prias.id, file = "prias_id.csv", row.names = F)


plot(survfit(Surv(year_discontinued, event_indicator)~1, conf.type="log-log", data=prias.id), mark.time = T)


#Ignore patient 957, no gleason for him
coxModel = coxph(Surv(year_discontinued, event_indicator) ~ poly(Age,2)[,1] + poly(Age,2)[,2], 
                 data=prias.id[!(prias.id$P_ID %in% c(957, 5102)),], x = T, model = T)
anova(coxModel)

#Cox snell residuals
coxsnellres=prias.id$event_indicator-resid(coxModel,type="martingale")
fitres=survfit(coxph(Surv(coxsnellres,prias.id$event_indicator)~1,method='breslow'),type='aalen')
plot(fitres$time,-log(fitres$surv),type='s',xlab='Cox-Snell Residuals', 
     ylab='Estimated Cumulative Hazard Function')
abline(0,1,col='red',lty=2)

