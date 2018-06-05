library(JMbayes)
library(splines)
library(ggplot2)

source("src/decision_analytic/simulationStudy/semiParamBH.R")

getTheoreticalHazard = function(jointModelData, timePoint){
  
  getBaselineHazard = function(weibullScale, weibullShape, times){
    return((weibullShape/weibullScale)*(times/weibullScale)^(weibullShape-1))
  }
  
  totalUnderLyingHazards = length(unique(jointModelData$trainingData$trainingDs.id$progression_speed))
  
  theoreticalHazard = sapply(1:totalUnderLyingHazards, function(index){
    progression_speed = unique(jointModelData$trainingData$trainingDs.id$progression_speed)[index]
    getBaselineHazard(weibullScale = jointModelData$weibullScales[progression_speed], 
                      weibullShape = jointModelData$weibullShapes[progression_speed], times = timePoint)
  })
  
  unscaledWeights = sapply(1:totalUnderLyingHazards, function(index){
    progression_speed = unique(jointModelData$trainingData$trainingDs.id$progression_speed)[index]
      weibullScale = jointModelData$weibullScales[progression_speed]
      weibullShape = jointModelData$weibullShapes[progression_speed]
      return(exp(-(timePoint/weibullScale)^weibullShape))
  })
  weights = unscaledWeights/sum(unscaledWeights)
  
  return(sum(theoreticalHazard * weights))
}

getSplineLogBaselineHazard = function(jointModelData, timePoint){
  splineDesign(jointModelData$mvJoint_psa_simDs$control$knots, timePoint, 
          ord = jointModelData$mvJoint_psa_simDs$control$ordSpline, outer.ok = T) %*% jointModelData$mvJoint_psa_simDs$statistics$postwMeans$Bs_gammas
}

getSemiParametricLogBaselineHazard = function(jointModelData, semiParamCoeffs, knots, timePoint){
  semiParamHazard = semiParamCoeffs[max(which(knots<=timePoint))]
}

times = seq(0, 10, length.out = 500)

savedFiles = list.files(path = "Rdata/decision_analytic/Simulation/Mixed/", full.names = T)
logHazardMatrix = matrix(ncol=length(times), nrow=length(savedFiles))

for(i in 1:length(savedFiles)){
  rm(jointModelData)
  load(savedFiles[i])
  print(paste("Reading the file number:", i))
  
  semiParamCoeffs = c(semiParamBH(jointModelData$survModel_simDs, length.knots = 20)[[1]]$coefficients[-c(1,2)], 0)
  knots = semiParamBH(jointModelData$survModel_simDs, length.knots = 20)[[2]]
  
  #logHazardMatrix[i, ] = sapply(times, getSemiParametricLogBaselineHazard, jointModelData=jointModelData, semiParamCoeffs=semiParamCoeffs, knots=knots)
  logHazardMatrix[i, ] = sapply(times, getSplineLogBaselineHazard, jointModelData=jointModelData)
    
  print(paste("******** End working on Data Set: ", i, "*******"))
}

mixtureTheoreticalLogHazard = log(sapply(times, getTheoreticalHazard, jointModelData=jointModelData))

fittedLogHazardMean = (apply((logHazardMatrix), 2, mean, na.rm=T))
fittedLogHazardLower = (apply((logHazardMatrix), 2, quantile, probs=c(0.025), na.rm=T))
fittedLogHazardUpper = (apply((logHazardMatrix), 2, quantile, probs=c(0.975), na.rm=T))

df = data.frame(logHazard = c(mixtureTheoreticalLogHazard, fittedLogHazardMean),
                "Baseline Hazard"=rep(c("Theoretical log (baseline hazard)", 
                                        "Fitted log(baseline hazard)"), each=500),
                confIntLow=rep(fittedLogHazardLower,2),
                confIntHigh=rep(fittedLogHazardUpper,2),
                time = rep(times,2))

ggplot(data=df[df$time > 2 & df$time<=10,]) + 
  geom_line(aes(x=time, y=logHazard, linetype=Baseline.Hazard)) + 
  geom_ribbon(aes(x=time, ymin=(confIntLow), ymax=(confIntHigh)), fill = "grey", alpha=0.4) + 
  scale_x_continuous(breaks = c(0.1, seq(1, 10, 1))) + 
  scale_linetype_manual(values=c("twodash", "solid")) +
  theme(legend.title=element_blank(), legend.position = "top",
        text = element_text(size=11), axis.text=element_text(size=11),
        legend.text = element_text(size=11)) + 
  xlab("Time (years)") + ylab("log (baseline hazard)")

ggsave(file="report/pers_schedule/biometrics_submission/images/sim_study/baseline_hazard.eps", 
       width=8.27*0.75, height=8.27*0.75, device=cairo_ps)
