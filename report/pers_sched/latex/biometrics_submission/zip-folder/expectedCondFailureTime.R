expectedCondFailureTime = function(object, newdata, idVar = "id", last.time=NULL, 
                                   maxPossibleFailureTime = NULL){
  
  if (!(inherits(object, "JMbayes")|inherits(object, "mvJMbayes")))
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