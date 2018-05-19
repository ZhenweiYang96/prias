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

source("src/pers_sched/replaceMCMCContents.R")
source("src/pers_sched/rocJM_mod.R")
source("../JMBayes/Anirudh/dev/grid_arrange_shared_legend.R")
source("../JMBayes/Anirudh/dev/expectedCondFailureTime.R")

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
#Either of the two models can be NULL
#use jointModelBayes object not mvJointModelBayes
plotLog2PSAJMFitT3VsNorm = function(jmfitNorm, jmfitT3, 
                                    patientId, maxRowNum=NA, 
                                    modelNames=c("Fitted (Normally distributed errors)","Fitted (t-distributed df=3, errors)")){
  
  ds = prias_long[prias_long$P_ID == patientId, ]
  
  if(!is.na(maxRowNum)){
    ds=ds[1:maxRowNum,]
  }
  
  lastBiopsyTime = tail(ds$visitTimeYears[!is.na(ds$gleason)],1)
  
  if(!is.null(jmfitNorm)){
    predNorm = predict(jmfitNorm, ds[!is.na(ds$psa),], type="Subject", idVar="P_ID", 
                       last.time = lastBiopsyTime,
                       FtTimes = seq(min(ds$visitTimeYears), max(ds$visitTimeYears), length.out = 50))
  }
  
  if(!is.null(jmfitT3)){
    predT3 = predict(jmfitT3, ds[!is.na(ds$psa),], type="Subject", idVar="P_ID", 
                     FtTimes = seq(min(ds$visitTimeYears), max(ds$visitTimeYears), length.out = 50))
  }
  
  ds = ds[,c("visitTimeYears", "log2psa")]
  ds$type="Observed"
  
  if(!is.null(jmfitNorm)){
    ds = rbind(ds, data.frame("visitTimeYears"=attributes(predNorm)$time.to.pred[[1]], 
                              "log2psa"=c(predNorm), type=modelNames[1]))
  }
  
  if(!is.null(jmfitT3)){
    ds = rbind(ds, data.frame("visitTimeYears"=attributes(predT3)$time.to.pred[[1]], 
                              "log2psa"=c(predT3), type=modelNames[2]))
  }
  
  breaks = if(round(max(ds$visitTimeYears))>6){
    seq(0, round(max(ds$visitTimeYears)), by=2)
  }else if(round(max(ds$visitTimeYears))>2){
    seq(0, round(max(ds$visitTimeYears)), by=1)
  }else if(round(max(ds$visitTimeYears))==2){
    seq(0, round(max(ds$visitTimeYears)), length.out=5)
  }else{
    seq(0, round(max(ds$visitTimeYears)), length.out=3)
  }
  
  linetypes = c("dotted", "dashed", "solid", "solid")
  
  if(is.null(jmfitNorm) | is.null(jmfitT3)){
    linetypes = c("dotted", "solid", "solid")  
  }
  
  p=ggplot(data=ds[!is.na(ds$log2psa),]) + 
    geom_line(data=ds[ds$type!="Observed",], aes(x=visitTimeYears, y=log2psa, linetype=type)) + 
    geom_point(size=4, data=ds[ds$type=="Observed",], aes(x=visitTimeYears, y=log2psa)) +
    scale_x_continuous(breaks=breaks) + 
    geom_vline(aes(xintercept = lastBiopsyTime, linetype="Latest biopsy"), show.legend=F) + 
    scale_linetype_manual(values=linetypes) + 
    theme(text = element_text(size=13), axis.text=element_text(size=13),
          legend.text = element_text(size=13), legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5, size=13)) +
    xlab("Time (years)") + ylab(expression('log'[2]*'(PSA)'))
  
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
        xlab("Time (Years)") + 
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
      xlab("Time (Years)") + 
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
        xlab("Time (Years)") + 
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
      xlab("Time (Years)") + 
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

plotDynamicSurvProb = function(pid, fittedJointModel, maxVisitTime, futureTimeDt = 3){
  #Do not use psa_data_set here. That was only created to have a data set of non NA PSA's
  
  patientDs = psa_data_set[psa_data_set$P_ID %in% pid & psa_data_set$visitTimeYears<=maxVisitTime,]
  
  lastPSATime = max(patientDs$visitTimeYears)
  lastBiopsyTime = getLastBiopsyTime(pid, upperLimitTime = maxVisitTime)
  #lastBiopsyTime = 0
  
  futureTimes = seq(lastBiopsyTime, lastBiopsyTime + futureTimeDt, 0.1)
  
  sfit = survfitJM(fittedJointModel, patientDs, idVar="P_ID", survTimes = futureTimes, last.time = lastBiopsyTime)
  
  longprof = predict(fittedJointModel, patientDs, type = "Subject",
                     interval = "confidence", return = TRUE, idVar="P_ID", FtTimes = futureTimes)
  
  longprof$log2psa[(nrow(patientDs)+1):nrow(longprof)] = NA
  longprof$survMean = rep(NA, nrow(longprof))
  longprof$survLow = rep(NA, nrow(longprof))
  longprof$survUp = rep(NA, nrow(longprof))
  
  ymin = min(c(longprof[longprof$visitTimeYears<=lastPSATime,]$pred, longprof[longprof$visitTimeYears<=lastPSATime,]$log2psa), na.rm = T)
  ymax = max(c(longprof[longprof$visitTimeYears<=lastPSATime,]$pred, longprof[longprof$visitTimeYears<=lastPSATime,]$log2psa), na.rm = T)
  
  #subsetting twice because there are two rows for the last time, and -1 to remove one of those two rows
  longprof[(nrow(patientDs)+1):nrow(longprof), c("survMean", "survLow", "survUp")] =
    (rbind(c(1,1,1),sfit$summaries[[1]][, c("Mean", "Lower", "Upper")]) * (ymax-ymin) + ymin)
  
  p=ggplot() +
    geom_line(data = longprof[longprof$visitTimeYears<=lastPSATime,], aes(x = visitTimeYears, y=pred), linetype="dashed") +
    geom_point(size=4, data = longprof[longprof$visitTimeYears<=lastPSATime,], aes(y=log2psa, x=visitTimeYears)) +
    geom_vline(xintercept = lastBiopsyTime, linetype="solid") +
    geom_line(data = longprof[(nrow(patientDs)+1):nrow(longprof),], aes(x=visitTimeYears, y=survMean)) +
    geom_ribbon(data = longprof[(nrow(patientDs)+1):nrow(longprof),], aes(ymin=survLow, ymax=survUp, x= visitTimeYears), fill="grey", alpha=0.5) +
    xlab("Time (years)") + ylab(expression('log'[2]*'(PSA)')) +
    theme(text = element_text(size=25), axis.text=element_text(size=25))  +
    scale_y_continuous(limits = c(ymin, ymax),sec.axis = sec_axis(~(.-ymin)/(ymax-ymin), name = "Dynamic survival probability"))
  
  return(p)
}

plotEwoutMeetingPlot = function(pid, fittedJointModel, maxVisitTime, futureTimeDt = 3){
  pp = plotDynamicSurvProb(pid, fittedJointModel, maxVisitTime, futureTimeDt)
  lastBiopsyTime = getLastBiopsyTime(pid, upperLimitTime = maxVisitTime)
  #lastBiopsyTime = 0
  
  patientDs = psa_data_set[psa_data_set$P_ID %in% pid & psa_data_set$visitTimeYears<=maxVisitTime,]
  lastPSATime = max(patientDs$visitTimeYears)
  
  futureTimes = seq(lastBiopsyTime, lastBiopsyTime + futureTimeDt, 0.1)
  sfit = survfitJM(fittedJointModel, patientDs, idVar="P_ID", survTimes = futureTimes, last.time = lastBiopsyTime)
  
  medianSurv = sfit$summaries[[1]][,"Median"]
  survTimes = sfit$summaries[[1]][,"times"]
  survTimeCutoff = c(survTimes[which(abs(medianSurv-0.9)==min(abs(medianSurv-0.9)))])
  
  #expecFailTime = expectedCondFailureTime(fittedJointModel, patientDs, 
                                          #idVar = "P_ID", last.time = lastBiopsyTime, maxPossibleFailureTime = 20)
  
  #print(expecFailTime)
  pp = pp + geom_text(aes(x=lastBiopsyTime, y=quantile(patientDs$log2psa, probs=0.25), label="Last Biopsy (Time = t)"), size=4, angle=90, nudge_x = -0.25)+
    geom_vline(xintercept = lastPSATime, linetype="solid") +
   geom_text(aes(x=lastPSATime, y=quantile(patientDs$log2psa, probs=0.25), label="Current Visit (Time = s)"), size=4, angle=90, nudge_x = -0.25) + 
    geom_vline(xintercept = survTimeCutoff, linetype="solid") +
    geom_text(aes(x=survTimeCutoff, y=quantile(patientDs$log2psa, probs=0.25), label="Survival Prob = 0.9"), size=4, angle=90, nudge_x = -0.25)
  
  print(pp)
}

#Chapter 5, section 5.2
#shiny::runGitHub("Repeated_Measurements", "drizopoulos")