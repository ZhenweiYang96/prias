library(survival)
library(survminer)
library(rms)

fit <- survfit(Surv(right_cens_time, reclassification) ~ 1, 
               data=prias_final.id)
ggsurvplot(fit, risk.table = TRUE)

survModel_rightCens = coxph(Surv(right_cens_time, reclassification) ~ age, 
                            data=prias_final.id, x = T, model = T)


fit=cph(Surv(right_cens_time, reclassification) ~ age,
        data=prias_final.id, surv=TRUE, x=T, y=T, u=1, time.inc = 1)

cal.KM <- calibrate.cph(fit, cmethod='KM',u=1)
# cal.hare=calibrate(fit,u=100,cmethod='hare',B=20)
# plot(cal.hare)
plot(cal.KM,add=TRUE)