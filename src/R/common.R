library(ggplot2)
library(ggmcmc)
library(coda)
library(parallel)
library(doParallel)
library(survival)
library(splines)
library(nlme)
library(JMbayes)

source("src/R/replaceMCMCContents.R")
source("src/R/rocJM_mod.R")
source("../JMBayes/Anirudh/dev/multiplot.R")

ticksX = function(from=0, max, by, labels=waiver(), extraLabels=NA){
  scale_x_continuous(breaks = c(seq(from, max, by = by), extraLabels), labels = labels)
}

ticksY = function(from=0, max, by, labels = waiver()){
  scale_y_continuous(breaks = round(seq(from, max, by = by),3), labels=waiver())
}

sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}

plotRandomProfile = function(count=1, fitted=F){
  pid_sample = sample(x = unique(training_psa_data_set$P_ID), size = count)
  print(pid_sample)
  plot<-ggplot(data=training_psa_data_set[training_psa_data_set$P_ID %in% pid_sample,], aes(x=visitTimeYears, y=log(psa, base = 2))) + 
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
      temp[,c(3,4,5)] = exp(temp[,c(3,4,5)])
    }
    temp
  })
  
  if(individually==T){
    return(lapply(plotData, function(data){
      plot = ggplot(data=data)
      if(transformPSA==F){
        plot = plot + geom_ribbon(aes(x=visitTimeYears, ymin=low, ymax=upp), fill = "grey", alpha=0.4)
      } 
      plot = plot + geom_line(aes(y=pred, x=visitTimeYears), color="black") +
        ticksX(from=0, max = 10, by = 1) + 
        ticksY(from=0, 25, by = if(transformPSA==F){0.1}else{0.5}) + 
        xlab("Follow-up time (Years)") + 
        ylab(expression('log'[2]*'(PSA)'))
      print(plot)
      return(plot)
    }))
  }else{
    newPlotData = do.call(rbind, plotData)
    newPlotData$model = rep(paste("Model",1:length(plotData)), each=nrow(newDF))
    
    plot = ggplot(data=newPlotData)
    if(transformPSA==F){
      plot = plot + geom_ribbon(aes(x=visitTimeYears, ymin=low, ymax=upp, group=model), 
                                fill = "grey", alpha=0.4)
    } 
    plot = plot + geom_line(aes(y=pred, x=visitTimeYears, group=model, color=model)) +
      ticksX(from=0, max = 10, by = 1) + 
      ticksY(from=0, 25, by = if(transformPSA==F){0.1}else{0.5}) + 
      xlab("Follow-up time (Years)") + 
      ylab(expression('log'[2]*'(PSA)'))
    print(plot)
    
    return(plot)
  }
}

fitUnivaritePSAModel = function(fixedSplineKnots=c(0.1,0.5, 4), randomSplineKnots=c(0.1), 
                                boundaryKnots=range(psa_data_set$visitTimeYears), method="ML"){

  fixedFormula = as.formula(paste("log2psa ~ I(Age - 70) + I((Age - 70)^2) + ",
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

getLastBiopsyTime = function(pid, upperLimitTime = Inf){
  temp = prias_long[prias_long$P_ID %in% pid & prias_long$visitTimeYears<=upperLimitTime,][, c("visitTimeYears", "gleason")]
  lastBiopsyTime = max(temp[complete.cases(temp),]$visitTimeYears)
  return(lastBiopsyTime)
}

plotDynamicSurvProb = function(pid, fittedJointModel, futureTimeDt = 3){
  #Do not use psa_data_set here. That was only created to have a data set of non NA PSA's
  lastBiopsyTime = getLastBiopsyTime(pid)
  
  patientDs = psa_data_set[psa_data_set$P_ID %in% pid,]
  lastPSATime = max(patientDs$visitTimeYears)
  futureTimes = seq(lastPSATime, lastPSATime + futureTimeDt, 0.1)
  
  sfit = survfitJM(fittedJointModel, patientDs, idVar="P_ID", survTimes = futureTimes, last.time = lastBiopsyTime)
  
  longprof = predict(fittedJointModel, patientDs, type = "Subject",
                     interval = "confidence", return = TRUE, idVar="P_ID", FtTimes = futureTimes)
  
  longprof[longprof$visitTimeYears>lastPSATime,]$logpsa1=NA
  longprof$survMean = rep(NA, nrow(longprof))
  longprof$survLow = rep(NA, nrow(longprof))
  longprof$survUp = rep(NA, nrow(longprof))
  
  ymin = min(c(longprof[longprof$visitTimeYears<=lastPSATime,]$pred, longprof[longprof$visitTimeYears<=lastPSATime,]$logpsa1))
  ymax = max(c(longprof[longprof$visitTimeYears<=lastPSATime,]$pred, longprof[longprof$visitTimeYears<=lastPSATime,]$logpsa1))
  
  #subsetting twice because there are two rows for the last time, and -1 to remove one of those two rows
  longprof[longprof$visitTimeYears>=lastPSATime, c("survMean", "survLow", "survUp")][-1,] =
    (sfit$summaries[[1]][, c("Mean", "Lower", "Upper")] * (ymax-ymin) + ymin)

  p=ggplot() +
    geom_line(data = longprof[longprof$visitTimeYears<=lastPSATime,], aes(x = visitTimeYears, y=pred)) +
    geom_point(data = longprof[longprof$visitTimeYears<=lastPSATime,], aes(y=logpsa1, x=visitTimeYears), colour="red", alpha=0.4) +
    geom_vline(xintercept = lastPSATime, linetype="dotted") +
    geom_line(data = longprof, aes(x=visitTimeYears, y=survMean)) +
    geom_ribbon(data = longprof, aes(ymin=survLow, ymax=survUp, x= visitTimeYears), fill="grey", alpha=0.5) +
    xlab("Time (years)") + ylab(expression('log'[2]*'(PSA)')) +
    scale_y_continuous(limits = c(ymin, ymax),breaks = round(seq(ymin, ymax, 0.25),2), 
                       sec.axis = sec_axis(~(.-ymin)/(ymax-ymin), name = "Dynamic survival probability"))

  print(p)
}


#Chapter 5, section 5.2
#shiny::runGitHub("Repeated_Measurements", "drizopoulos")