library(ggplot2)
library(ggmcmc)
library(coda)
library(parallel)
library(doParallel)
library(survival)
library(splines)
library(nlme)
library(JMbayes)

ticksX = function(from=0, max, by, labels=waiver()){
  scale_x_continuous(breaks = seq(from, max, by = by), labels = labels)
}

ticksY = function(from=0, max, by, labels = waiver()){
  scale_y_continuous(breaks = seq(from, max, by = by), labels=waiver())
}

sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}

plotRandomProfile = function(count=1, fitted=F){
  pid_sample = sample(x = unique(psa_data_set$P_ID), size = count)
  plot<-ggplot(data=psa_data_set[psa_data_set$P_ID %in% pid_sample,], aes(x=visitTimeYears, y=log(psa + 1, base = 2))) + 
    geom_line(aes(group=P_ID))
  if(fitted==T){
    plot + geom_line(aes(y=fitted, x=visitTimeYears, color=P_ID, group=P_ID)) 
  }else{
    plot
  }
}

#######################################################################
# the following function creates the predicted values
# and the 95% CIs
#######################################################################
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

plotPSAFittedCurve = function(models, transformPSA=F, individually=T){
  newDF <- with(psa_data_set, expand.grid(visitTimeYears = seq(0, 10, by = 0.5),
                                          Age = 70))
  plotData = lapply(models, function(model){
    temp = effectPlotData(model, newDF, psa_data_set)
    
    if(transformPSA==T){
      temp[,c(3,4,5)] = exp(temp[,c(3,4,5)]) - 1
    }
    temp
  })
  
  if(individually==T){
    lapply(plotData, function(data){
      plot = ggplot(data=data)
      if(transformPSA==F){
        plot = plot + geom_ribbon(aes(x=visitTimeYears, ymin=low, ymax=upp), fill = "grey", alpha=0.4)
      } 
      plot = plot + geom_line(aes(y=pred, x=visitTimeYears), color="black") +
        ticksX(from=0, max = 10, by = 1, labels = c(0:10)) + 
        ticksY(from=0, 25, by = if(transformPSA==F){0.1}else{0.5}) + 
        xlab("Follow-up time (Years)") + 
        ylab("PSA")
      print(plot)
    })
  }else{
    newPlotData = do.call(rbind, plotData)
    newPlotData$model = rep(paste("Model",1:length(plotData)), each=nrow(newDF))
    
    plot = ggplot(data=newPlotData)
    if(transformPSA==F){
      plot = plot + geom_ribbon(aes(x=visitTimeYears, ymin=low, ymax=upp, group=model), 
                                fill = "grey", alpha=0.4)
    } 
    plot = plot + geom_line(aes(y=pred, x=visitTimeYears, group=model, color=model)) +
      ticksX(from=0, max = 10, by = 1, labels = c(0:10)) + 
      ticksY(from=0, 25, by = if(transformPSA==F){0.1}else{0.5}) + 
      xlab("Follow-up time (Years)") + 
      ylab("PSA")
    print(plot)
  }
}

pltoDynamicPredictions = function(pid, fittedJointModel, survival=T, longitudinal=T){
    ND = psa_data_set[psa_data_set$P_ID %in% pid,]
    futureTimes = seq(max(ND$visitTimeYears), min(10,(max(ND$visitTimeYears) + 3)), 0.1)
    
    sfit.patient_temp = survfitJM(fittedJointModel, ND, idVar="P_ID", survTimes = futureTimes)
    plot(sfit.patient_temp, estimator="mean", include.y=T, conf.int=T, fill.area=T, col.area="lightgrey", main=paste("P_ID =", pid), xlab="Time (years)")
    
    longprof = predict(fittedJointModel, ND, type = "Subject",
                       interval = "confidence", return = TRUE, idVar="P_ID", FtTimes = futureTimes)
    last.time <- with(longprof, visitTimeYears[!is.na(low)][1])
    lattice::xyplot(pred + low + upp ~ visitTimeYears, data = longprof, type = "l",
                    lty = c(1, 2, 2), col = c(2, 1, 1), abline = list(v = last.time, lty = 3),
                    xlab = "Time (years)", ylab = "Predicted log(PSA + 1)", main=paste("P_ID =",id))
}

fitUnivaritePSAModel = function(fixedSplineKnots=c(1,2,4), randomSplineKnots=c(1), 
                                boundaryKnots=range(psa_data_set$visitTimeYears), method="ML"){

  fixedFormula = as.formula(paste("logpsa1 ~ I(Age - 70) + I((Age - 70)^2) + ",
                  "ns(visitTimeYears, knots=c(", paste(fixedSplineKnots, collapse=", "), 
                  "), Boundary.knots=c(", paste(boundaryKnots, collapse=", "),"))", sep=""))
  
  randomFormula = as.formula(paste("~ns(visitTimeYears, knots=c(", paste(randomSplineKnots, collapse=", "), 
                        "), Boundary.knots=c(", paste(boundaryKnots, collapse=", "),"))|P_ID", sep=""))
  
  model = lme(fixed=fixedFormula, random = randomFormula, 
      data=psa_data_set,
      control = lmeControl(opt = "optim", optimMethod = "L-BFGS-B"), 
      method = method)
  
  model$call$fixed = fixedFormula
  model$call$random = randomFormula
  
  print(anova(model, type="marginal"))
  plot = qplot(x=model$fitted[,2], y=model$residuals[,2], xlab="Fitted", 
        ylab="Residuals", geom=c("point", "smooth"))
  print(plot)
  
  return(model)
}


#Chapter 5, section 5.2
#shiny::runGitHub("Repeated_Measurements", "drizopoulos")