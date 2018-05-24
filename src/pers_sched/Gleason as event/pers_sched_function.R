library(JMbayes)
library(splines)

modifyScheduledBiopsyTime = function(proposedTime, curVisitTime, lastBiopsyTime){
  if(proposedTime < curVisitTime){
    if(curVisitTime - lastBiopsyTime <= 1){
      return(lastBiopsyTime + 1)
    }else{
      return(curVisitTime)
    }
  }else{
    if(proposedTime - lastBiopsyTime <= 1){
      return(lastBiopsyTime + 1)
    }else{
      return(proposedTime)
    }
  }
}

expectedCondFailureTime = function(object, newdata, idVar = "id", last.time=NULL, 
                                   maxPossibleFailureTime = NULL){
  
  if (!inherits(object, "JMbayes"))
    stop("Use only with 'JMbayes' objects.\n")
  if (!is.data.frame(newdata) || nrow(newdata) == 0L)
    stop("'newdata' must be a data.frame with more than one rows.\n")
  if (is.null(newdata[[idVar]]))
    stop("'idVar' not in 'newdata.\n'")
  
  dynamicPredProb = function(futureTimes, object, newdata, last.time, idVar){
    return(survfitJM(object, newdata, last.time = last.time,
                     idVar=idVar, survTimes = futureTimes)$summaries[[1]][, "Median"])
  }
  
  if(is.null(last.time)){
    last.time =  max(newdata[[object$timeVar]])
  }
  
  if(is.null(maxPossibleFailureTime)){
    maxPossibleFailureTime =  max(object$y$Time) * 1.5
  }
  
  last.time + integrate(dynamicPredProb, lower=last.time, 
                        upper=maxPossibleFailureTime, object=object, newdata=newdata, last.time = last.time,
                        idVar=idVar, rel.tol = 0.05)$value
}

nextBiopsyTime = function(object, newdata, idVar = "id", 
                          maxPossibleFailureTime = NULL, 
                          methodName="expectedFailureTime",
                          dynSurvCutoff = NA){

  if(is.null(maxPossibleFailureTime)){
    maxPossibleFailureTime =  max(object$y$Time) * 1.5
  }
  
  curVisitTime = tail(newdata$visitTimeYears, 1)
  lastBiopsyTime = tail(newdata$visitTimeYears[!is.na(newdata$gleason)],1)

  proposedBiopsyTime = NA
  if(methodName=="expectedFailureTime"){
    proposedBiopsyTime = expectedCondFailureTime(object, newdata, idVar, 
                                                 lastBiopsyTime, maxPossibleFailureTime)
  }else if(methodName %in% c("f1score", "youden", "fixedCutoff")){
    survTimes = seq(lastBiopsyTime, maxPossibleFailureTime, 0.1)
    survProbs = c(1,survfitJM(object, newdata[!is.na(newdata$psa),],
                              idVar=idVar, last.time = lastBiopsyTime,
                              survTimes = survTimes)$summaries[[1]][, "Median"])
    

    survProbCutoff = dynSurvCutoff
    if(methodName %in% c("f1score", "youden")){
      nearest_time_index = which(abs(dynamicCutoffTimes_PRIAS-lastBiopsyTime)==min(abs(dynamicCutoffTimes_PRIAS-lastBiopsyTime)))[1]
      survProbCutoff = cutoffValues_PRIAS[[nearest_time_index]][methodName]
    }
    
    if(is.na(survProbCutoff)){
      return(NA)
    }else{
      proposedBiopsyTime = survTimes[which(abs(survProbs-survProbCutoff)==min(abs(survProbs-survProbCutoff)))[1]]
    }
  }
  
  return(modifyScheduledBiopsyTime(proposedBiopsyTime, curVisitTime, lastBiopsyTime))
}