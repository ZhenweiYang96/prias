library(nlme)
library(splines)
library(ggplot2)

source("common.R")

psa_data_set =  prias_long[!is.na(prias_long$psa),]
psa_data_set = prias_long_survival_eligible[!is.na(prias_long_survival_eligible$psa),]

psa_data_set$P_ID = droplevels(psa_data_set$P_ID)
psa_data_set$logpsa1 = log(psa_data_set$psa + 1)
save.image(file = "psa_gl.Rdata")

plotRandomProfile = function(count=1, fitted=F){
  pid_sample = sample(x = unique(psa_data_set$P_ID), size = count)
  plot<-ggplot(data=psa_data_set[psa_data_set$P_ID %in% pid_sample,], aes(x=visitTimeYears, y=log(psa + 1))) + 
    geom_line(aes(group=P_ID))
  if(fitted==T){
    plot + geom_line(aes(y=fitted, x=visitTimeYears, color=P_ID, group=P_ID)) 
  }else{
    plot
  }
}

psa_data_set$residuals = rep(NA, nrow(psa_data_set))
psa_data_set$residuals[!is.na(psa_data_set$psa)] = residuals(lm(log(psa + 1)~poly(Age,2), data = psa_data_set))

pidList = unique(psa_data_set$P_ID)
randomSample = sample(x = pidList, replace = F, size = 150)
ggplot(data=psa_data_set[psa_data_set$P_ID %in% randomSample,], aes(x=visitTimeYears, y=log(psa+1))) + 
  geom_line(aes(group=P_ID)) + stat_smooth() + 
  ticksX(from=0,max(psa_data_set$visitTimeYears), 0.25) + xlab("Time (years)") + 
  ylab("log(psa + 1): Adjusted for Age")

#Since every box plot has 1 entry per person, 
#if the first few box plots have too many time points it means that it is not because
#people came frequently, it could also be due to the fact that people started dropping out later
#and so initially we see more concentration of time points.
#Per person the concentration is not more in the beginning, but overall marginally
#there is more concentration in the beginning
ggplot(data=prias_long[prias_long$P_ID,], aes(factor(visit_number),visitTimeYears)) + 
  geom_boxplot() + stat_summary(fun.data = function(x){
    return(c(y = 1.2, label = length(x))) 
  }, geom = "text", fun.y = median) + ticksY(max(prias_long$visitTimeYears), 0.05)

quantiles_time_years = quantile(prias_long$visitTimeYears, probs = c(0.25, 0.5, 0.75))
# 25%        50%        75% 
# 0.04849315 0.12136986 0.24849315 

#I would take equidistant time points for cutoff
prias_long = cbind(prias_long, polyage=poly(prias_long$Age,3))
prias_long = cbind(prias_long, nsFixed=ns(prias_long$visitTimeYears, knots=quantiles_time_years),
                   nsRandom=ns(prias_long$visitTimeYears, knots=quantiles_time_years[1]))
write.csv(prias_long, "prias_long.csv", row.names = F)

psaModel_big = lme(fixed=log(psa + 1)~poly(Age,2)[,1] + poly(Age,2)[,2] +
                 ns(visitTimeYears, knots=c(0.1, 0.35, 0.6, 0.85)),
                     random = ~ns(visitTimeYears, knots=c(0.1))|P_ID, 
                     data=prias_long[!is.na(prias_long$psa),],
               control = lmeControl(opt = "optim", optimMethod = "L-BFGS-B"), 
               method = "ML")

anova(psaModel_big, type="marginal")

psaModel_1 = lme(fixed=log(psa + 1)~poly(Age,2)[,1] + poly(Age,2)[,2] +
                     ns(visitTimeYears, knots=c(0.1, 0.35, 0.6, 0.85)) * poly(Age,2)[,1],
                   random = ~ns(visitTimeYears, knots=c(0.1))|P_ID, 
                   data=prias_long[!is.na(prias_long$psa),],
                   control = lmeControl(opt = "optim", optimMethod = "L-BFGS-B"), 
                   method = "ML")

prias_long$fitted[!is.na(prias_long$psa)] = psaModel_1$fitted[,2]
prias_long$residuals[!is.na(prias_long$psa)] = psaModel_1$residuals[,2]

plotRandomProfile(1,T)

psaModel_2 = lme(fixed=log(psa + 1)~poly(Age,2)[,1] + poly(Age,2)[,2] + 
                   ns(visitTimeYears, knots=c(0.03, 0.09, 0.16, 0.30))*poly(Age,2)[,1],
                 random = ~ns(visitTimeYears, knots=c(0.1))|P_ID, 
                 data=prias_long[!is.na(prias_long$psa),],
                 control = lmeControl(opt = "optim", optimMethod = "L-BFGS-B"), 
                 method = "ML")

anova(psaModel_1, psaModel_2)

psaModel_3_1 = lme(fixed=log(psa + 1)~poly(Age,2)[,1] + poly(Age,2)[,2] + 
                   ns(visitTimeYears, knots=c(0.1, 0.2, 0.4))*poly(Age,2)[,1],
                 random = ~ns(visitTimeYears, knots=c(0.1))|P_ID, 
                 data=prias_long[!is.na(prias_long$psa) & !(prias_long$P_ID %in% c(957, 5102)),],
                 control = lmeControl(opt = "optim", optimMethod = "L-BFGS-B"), 
                 method = "ML")

prias_long$fitted[!is.na(prias_long$psa)] = psaModel_3_1$fitted[,2]
prias_long$residuals[!is.na(prias_long$psa)] = psaModel_3_1$residuals[,2]

anova(psaModel_1, psaModel_2, psaModel_3_1)

psaModel_4 = lme(fixed=log(psa + 1)~poly(Age,2)[,1] + poly(Age,2)[,2] + 
                   ns(visitTimeYears, knots=c(0.03, 0.09, 0.16, 0.30, 0.7))*poly(Age,2)[,1],
                 random = ~ns(visitTimeYears, knots=c(0.03))|P_ID, 
                 data=prias_long[!is.na(prias_long$psa),],
                 control = lmeControl(opt = "optim", optimMethod = "L-BFGS-B"), 
                 method = "ML")

anova(psaModel_1, psaModel_2, psaModel_4, psaModel_3 )


psaModel_5 = lme(fixed=log(psa + 1)~poly(Age,2)[,1] + poly(Age,2)[,2] + 
                   ns(visitTimeYears, knots=c(0.03, 0.09, 0.16, 0.30, 0.7))*poly(Age,2)[,1],
                 random = ~ns(visitTimeYears, knots=c(0.03, 0.09))|P_ID, 
                 data=prias_long[!is.na(prias_long$psa),],
                 control = lmeControl(opt = "optim", optimMethod = "L-BFGS-B"), 
                 method = "ML")

anova(psaModel_1, psaModel_2, psaModel_4, psaModel_3 )

###########################
###########################
psaModel_linear = lme(fixed=logpsa1~I(Age - 70) +  I((Age - 70)^2) + 
                        visitTimeYears,
                     random = ~visitTimeYears|P_ID, 
                     data=psa_data_set,
                     control = lmeControl(opt = "optim", optimMethod = "L-BFGS-B"))

psaModel_spline = lme(fixed=logpsa1~I(Age - 70) +  I((Age - 70)^2) + 
                   ns(visitTimeYears, knots=c(1, 2, 4)),
                 random = ~ns(visitTimeYears, knots=c(0.1))|P_ID, 
                 data=psa_data_set,
                 control = lmeControl(opt = "optim", optimMethod = "L-BFGS-B"))

mvglm_psa_linear = mvglmer(list(
  logpsa1 ~  I(Age - 70) +  I((Age - 70)^2) + visitTimeYears + 
    (visitTimeYears|P_ID)), 
  data =psa_data_set, families = list(gaussian))

mvJoint_psa_linear = mvJointModelBayes(mvglm_psa_linear, coxModel, timeVar = "visitTimeYears") 
save.image("psa_gl.Rdata")                                   

mvglm_psa = mvglmer(list(
  logpsa1 ~  I(Age - 70) +  I((Age - 70)^2) + 
    ns(visitTimeYears, knots=c(1, 2, 4), Boundary.knots=c(0,15)) + 
    (ns(visitTimeYears, knots=c(1), Boundary.knots=c(0,15))|P_ID)), 
  data =psa_data_set, families = list(gaussian))

save.image("psa_gl.Rdata")

mvJoint_psa = mvJointModelBayes(mvglm_psa, coxModel, timeVar = "visitTimeYears", 
                                   Formulas = list("logpsa1" = "value",
                                                   "logpsa1" = list(fixed = ~ 0 + dns(visitTimeYears, knots=c(1, 2, 4), Boundary.knots=c(0,15)),
                                                                    random=~0 + dns(visitTimeYears, knots=c(1), Boundary.knots=c(0,15)),
                                                                    indFixed = 4:7, indRandom=2:3, name = "slope")))
save.image("psa_gl.Rdata")

joint_psa = jointModelBayes(psaModel_spline, coxModel, timeVar = "visitTimeYears",
                                         param = "td-both",
                                         extraForm = list(fixed = ~ 0 + dns(visitTimeYears, knots=c(1, 2, 4), Boundary.knots=c(0,15)),
                                                          random=~0 + dns(visitTimeYears, knots=c(1), Boundary.knots=c(0,15)),
                                                          indFixed = 4:7, indRandom=2:3))

joint_psa_replaced = replaceMCMCContents(fromObj = mvJoint_psa, toObj = joint_psa)
save.image("psa_gl.Rdata")

#######################################################################
# the following function creates the predicted values
# and the 95% CIs
effectPlotData <- function (object, newdata, orig_data) {
  form <- formula(object)
  namesVars <- all.vars(form)
  betas <- if (!inherits(object, "lme")) coef(object) else fixef(object)
  V <- if (inherits(object, "geeglm")) object$geese$vbeta else vcov(object)
  orig_data <- orig_data[complete.cases(orig_data[namesVars]), ]
  Terms <- delete.response(terms(form))
  mfX <- model.frame(Terms, data = orig_data)
  Terms_new <- attr(mfX, "terms")
  mfX_new <- model.frame(Terms_new, newdata, xlev = .getXlevels(Terms, mfX))
  X <- model.matrix(Terms_new, mfX_new)
  pred <- c(X %*% betas)
  ses <- sqrt(diag(X %*% V %*% t(X)))
  newdata$pred <- pred
  newdata$low <- pred - 1.96 * ses
  newdata$upp <- pred + 1.96 * ses
  newdata
}

# the data frame that contains the combination of values to
# create the plot
newDF <- with(psa_data_set, expand.grid(visitTimeYears = seq(0, 10, by = 0.5),
                                 Age = 70))
plotData = effectPlotData(psaModel_spline, newDF, psa_data_set)
  
plotData[,c(3,4,5)] = exp(plotData[,c(3,4,5)]) - 1

#effects plot
ggplot(data=plotData) + geom_line(aes(y=pred, x=visitTimeYears)) + 
  ticksX(from=0, max = 10, by = 1, labels = c(0:10)) + ticksY(from=0, 10, by = 0.1) +xlab("Follow-up time (Years)") + 
ylab("PSA")

#Hazard ratio plot
sim_psa_slope = seq(1, 1 * 5, by=1)
simhr = exp(4.583786593 * (sim_psa_slope-sim_psa_slope[1]))

qplot(x = sim_psa_slope, y=simhr, geom = "line", ylab = "Hazard Ratio", xlab = "PSA Multiplier") + 
  ticksX(from=sim_psa_slope[1], max = max(sim_psa_slope), by = sim_psa_slope[1]*2) + ticksY(from=0, 4.5, by = 0.25)

#survival plot: Dynamic predictions
idList = unique(droplevels(prias.id[prias.id$progressed==0 & prias.id$DiscontinuedYesNo==0,]$P_ID))

idList = c(1817, 332, 1370)
for(id in idList){
  id = sample(idList, 1)
  id = idList[4]
  ND = psa_data_set[psa_data_set$P_ID %in% id,]
  futureTimes = seq(max(ND$visitTimeYears), min(12,(max(ND$visitTimeYears) + 3)), 0.1)
  
  sfit.patient_temp = survfitJM(joint_psa_replaced, ND, idVar="P_ID", survTimes = futureTimes)
  plot(sfit.patient_temp, estimator="mean", include.y=T, conf.int=T, fill.area=T, col.area="lightgrey", main=paste("P_ID =",id), xlab="Time (years)")
   
  longprof = predict(joint_psa_replaced, ND, type = "Subject",
                      interval = "confidence", return = TRUE, idVar="P_ID", FtTimes = futureTimes)
  last.time <- with(longprof, visitTimeYears[!is.na(low)][1])
   lattice::xyplot(pred + low + upp ~ visitTimeYears, data = longprof, type = "l",
                   lty = c(1, 2, 2), col = c(2, 1, 1), abline = list(v = last.time, lty = 3),
                   xlab = "Time (years)", ylab = "Predicted log(PSA + 1)", main=paste("P_ID =",id))
}