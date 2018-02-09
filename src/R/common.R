library(ggplot2)
library(ggmcmc)
library(coda)
library(parallel)
library(doParallel)
library(survival)
library(splines)
library(nlme)
library(JMbayes)
library(latex2exp)

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

plotProfile = function(fitted=F, patientRowNum=1){
  pid_sample = training.prias.id$P_ID[sample(1:length(training.prias.id), size = 1)]
  print(pid_sample)
  plot<-ggplot(data=training_psa_data_set[training_psa_data_set$P_ID %in% pid_sample,], aes(x=visitTimeYears, y=log(psa, base = 2))) + 
    geom_line(aes(group=P_ID))
  if(fitted==T){
    plot + geom_line(aes(y=fitted, x=visitTimeYears, color="new")) + 
      geom_line(aes(y=fitted_2, x=visitTimeYears, color="old")) 
  }else{
    plot
  }
}

#######################################################################
# the following function creates the predicted values
# and the 95% CIs
#######################################################################
effectPlotLmeData <- function (object, newdata, orig_data) {
  form <- formula(object)
  namesVars <- all.vars(form)
  
  betas <- if (!inherits(object, "lme")) coef(object) else fixef(object)
  
  #you can replace the contents of lme object by these coefficients
  #betas <- joint_psa_replaced_prias_t3$postMeans$betas
  #betas <- mvJoint_psa_spline_pt1pt54_pt1_tdboth_light$statistics$postMeans$betas1
  #betas <- mvJoint_psa_spline_pt1pt54_pt1_tdboth_t3$statistics$postMeans$betas1
  
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

effectPlotJMData <- function (mvJMbayesObject, lmeObject, newdata, orig_data) {
  form <- formula(lmeObject)
  namesVars <- all.vars(form)
  
  betas <- mvJMbayesObject$mcmc$betas1
  
  orig_data <- orig_data[complete.cases(orig_data[namesVars]), ]
  Terms <- delete.response(terms(form))
  mfX <- model.frame(Terms, data = orig_data)
  Terms_new <- attr(mfX, "terms")
  mfX_new <- model.frame(Terms_new, newdata, xlev = .getXlevels(Terms, mfX))
  X <- model.matrix(Terms_new, mfX_new)
  
  predMcmc <- betas %*% t(X)
  
  newdata$pred <- apply(predMcmc, 2, mean)
  newdata$low <- apply(predMcmc, 2, quantile, probs=0.025)
  newdata$upp <- apply(predMcmc, 2, quantile, probs=0.975)
  newdata
}

#not mvjoint model bayes
plotLog2PSAJMFit = function(jmfit, patientId){
  
  ds = pria[psa_data_set$P_ID == patientId, ]
  
  pred = predict(jmfit, ds, type="Subject", idVar="P_ID", 
                 FtTimes = seq(min(ds$visitTimeYears), max(ds$visitTimeYears), length.out = 50))
  
  ds = ds[,c("visitTimeYears", "log2psa")]
  ds$type="Observed"
  
  ds = rbind(ds, data.frame("visitTimeYears"=attributes(pred)$time.to.pred[[1]], 
                            "log2psa"=c(pred), type="Fitted"))
 
  
  breaks = if(round(max(ds$visitTimeYears))>2){
    seq(0, round(max(ds$visitTimeYears)), by=1)
  }else if(round(max(ds$visitTimeYears))==2){
    seq(0, round(max(ds$visitTimeYears)), length.out=5)
  }else{
    seq(0, round(max(ds$visitTimeYears)), length.out=3)
  }
    
  p=ggplot(data=ds) + 
    geom_line(aes(x=visitTimeYears, y=log2psa, linetype=type)) + 
    scale_linetype_manual(values=c("twodash", "solid")) + 
    scale_x_continuous(breaks=breaks) + 
    theme(text = element_text(size=11), axis.text=element_text(size=11),
          legend.text = element_text(size=11), legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5, size=13)) +
    xlab("Time(years)") + ylab(expression('log'[2]*'(PSA)'))
  
  print(p)
  
  return(p)
}


#Example usage
# plotMarginalPSAFittedCurveJM(list(mvJoint_psa_spline_pt1pt54_pt1_tdboth_light, 
# mvJoint_psa_spline_pt1pt54_pt1_tdboth_t3),
# lme_psa_spline_pt1pt54_pt1, individually = F, 
# modelNames = c("Normally distributed errors", "t-distributed (df=3) errors"))
plotMarginalPSAFittedCurveJM = function(models, lmeObject, transformPSA=F, individually=T, modelNames=NULL){
  newDF <- with(training_psa_data_set, expand.grid(visitTimeYears = seq(0, 10, by = 0.5),
                                                   Age = 70))
  plotData = lapply(models, function(model){
    temp = effectPlotJMData(model, lmeObject, newDF, training_psa_data_set)
    
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
        ylab(expression('Fitted marginal '*'log'[2]*'(PSA)')) + 
        theme(text = element_text(size=13), axis.text=element_text(size=13))
      print(plot)
      return(plot)
    }))
  }else{
    newPlotData = do.call(rbind, plotData)
    
    if(!is.null(modelNames)){
      newPlotData$Model = rep(modelNames, each=nrow(newDF))
    }else{
      newPlotData$Model = rep(paste("Model",1:length(plotData)), each=nrow(newDF))  
    }
    
    plot = ggplot(data=newPlotData)
    if(transformPSA==F){
      plot = plot + geom_ribbon(aes(x=visitTimeYears, ymin=low, ymax=upp, group=Model), 
                                fill = "grey", alpha=0.4)
    } 
    plot = plot + geom_line(aes(y=pred, x=visitTimeYears, group=Model, linetype=Model)) +
      ticksX(from=0, max = 10, by = 1) + 
      ticksY(from=0, 25, by = if(transformPSA==F){0.1}else{0.5}) + 
      xlab("Follow-up time (Years)") + 
      ylab(expression('Fitted marginal '*'log'[2]*'(PSA)')) + 
      theme(text = element_text(size=11), axis.text=element_text(size=11), 
            legend.position = "top", legend.title = element_blank(), 
            legend.text = element_text(size=11)) + 
      scale_linetype_manual(values=c("twodash", "solid")) 
    print(plot)
    
    return(plot)
  }
}

plotPSAFittedCurveLME = function(models, transformPSA=F, individually=T, modelNames=NULL){
  newDF <- with(training_psa_data_set, expand.grid(visitTimeYears = seq(0, 10, by = 0.5),
                                          Age = 70))
  plotData = lapply(models, function(model){
    temp = effectPlotData(model, newDF, training_psa_data_set)
    
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
        ylab(expression('log'[2]*'(PSA)')) + theme(text = element_text(size=13), axis.text=element_text(size=13))
      print(plot)
      return(plot)
    }))
  }else{
    newPlotData = do.call(rbind, plotData)
    
    if(!is.null(modelNames)){
      newPlotData$Model = rep(modelNames, each=nrow(newDF))
    }else{
      newPlotData$Model = rep(paste("Model",1:length(plotData)), each=nrow(newDF))  
    }
    
    plot = ggplot(data=newPlotData)
    if(transformPSA==F){
      plot = plot + geom_ribbon(aes(x=visitTimeYears, ymin=low, ymax=upp, group=Model), 
                                fill = "grey", alpha=0.4)
    } 
    plot = plot + geom_line(aes(y=pred, x=visitTimeYears, group=Model, linetype=Model)) +
      ticksX(from=0, max = 10, by = 1) + 
      ticksY(from=0, 25, by = if(transformPSA==F){0.1}else{0.5}) + 
      xlab("Follow-up time (Years)") + 
      ylab(expression('log'[2]*'(PSA)')) + 
      theme(text = element_text(size=11), axis.text=element_text(size=11), 
            legend.position = "top", legend.title = element_blank()) + 
      scale_linetype_manual(values=c("twodash", "solid")) 
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

getLastBiopsyTime = function(pid, lastnumber=1, upperLimitTime = Inf){
  temp = prias_long[prias_long$P_ID %in% pid & prias_long$visitTimeYears<=upperLimitTime,][, c("visitTimeYears", "gleason")]
  lastBiopsyTime = tail(temp[complete.cases(temp),]$visitTimeYears, lastnumber)[1]
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