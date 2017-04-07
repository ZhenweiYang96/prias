source("src/R/Simulation Study/ExpectedCondFailureTime.R")

plotTrueSurvival = function(patientId){
  
  time = 1:15
  survProb = sapply(time, function(t){
    survivalFunc(t, patientId)
  })
  
  byYAxis = (max(survProb) - min(survProb))/10
  
  qplot(x=time, y = survProb, geom = "line", xlab = "Time (years)", ylab = "Probability") + 
    ticksX(from=0, max = 15, by=1) + ticksY(min(survProb), max(survProb), by=byYAxis) + 
    ggtitle(patientId)
}

plotTrueLongitudinal = function(patientId){
  ggplot(data=simDs[simDs$P_ID == patientId, ], aes(x=visitTimeYears, y=logpsa1)) + 
    geom_line() + geom_point(color="red") + ggtitle(patientId) + xlab("Time (years)") + 
    ylab("log(PSA + 1)")
}

plotDynamicSurvival = function(patientId){
  ggplot(data=simTestDs[simTestDs$P_ID==patientId,]) + 
    geom_line(aes(x=visitTimeYears, y=fixed_pt5yr_survprob)) + xlab("Time (years)") + 
    ylab("Probability") + ggtitle(patientId)
}

simJointModel_replaced = replaceMCMCContents(mvJoint_psa_tdboth_training, jmbayes_psa_tdboth_training)

invDynSurvival <- function (t, u, patientDs) {
  u - survfitJM(simJointModel_replaced, patientDs, idVar="P_ID", survTimes = t)$summaries[[1]][1, "Mean"]
}

pDynSurvTime = function(survProb, patientDs){
  #Return the time at which the dynamic survival probability is say 90%
  
  Low = max(patientDs$visitTimeYears) + 1e-05
  Up <- 15
  tries  = 0
  
  repeat{
    tries = tries + 1
    Root <- try(uniroot(invDynSurvival, interval = c(Low, Up), 
                        u = survProb, patientDs = patientDs)$root, TRUE)
    
    if(inherits(Root, "try-error")){
      if(tries >= 5){
        return(NA)
      }else{
        Up = Up + 0.5    
      }
    }else{
      return(Root)
    }
  }
}

rLogPSA1 =  function(patientId, time){
  df_s = data.frame(visitTimeYears = time, Age = simDs.id$Age[patientId])
  
  xi_s_val = model.matrix(fixedValueFormula, df_s)
  zi_s_val = model.matrix(randomValueFormula, df_s)
  zib_val = zi_s_val %*% b[patientId, ]
  xBetaZb_s_value = xi_s_val %*% betas + zib_val
  
  sapply(xBetaZb_s_value, rnorm, n=1, sigma.y)
}