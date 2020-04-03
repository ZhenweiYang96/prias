library(JMbayes)
library(splines)
library(ggplot2)

semiParamBH = function (coxObject, knots = NULL, length.knots = 6) {
  Time <- coxObject$y[, 1]
  d <- coxObject$y[, 2]
  n <- length(Time)
  if (is.null(knots)) {
    Q <- length.knots + 1
    knots <- unique(quantile(Time, seq(0, 1, len = Q + 1), 
                             names = FALSE)[-c(1, Q + 1)])
    knots <- knots + 1e-06
    if (max(knots) > max(Time)) 
      knots[which.max(knots)] <- max(Time) - 1e-06
  }
  knots <- c(0, sort(knots), max(Time) + 1)
  Q <- length(knots) - 1
  ind <- findInterval(Time, knots, rightmost.closed = TRUE)
  D <- matrix(0, n, Q)
  D[cbind(seq_along(ind), ind)] <- 1
  D <- c(D * d)
  Tiq <- outer(Time, knots, pmin)
  T <- c(Tiq[, 2:(Q + 1)] - Tiq[, 1:Q])
  X <- coxObject$x[rep(seq_len(n), Q), ]
  ND <- suppressWarnings(data.frame(Time = T, D = D, X, xi = gl(Q, 
                                                                n), check.names = FALSE)[T > 0, ])
  return(list(glm(D ~ . + offset(log(Time)) - Time - 1, family = poisson, 
                  data = ND), knots))
}

getTheoreticalHazard = function(timePoint, jointModelData){
  
  weibullScales = jointModelData$weibullScales
  weibullShapes = jointModelData$weibullShapes
  progression_speeds = jointModelData$progression_speeds
  
  getBaselineHazard = function(weibullScale, weibullShape, times){
    return((weibullShape/weibullScale)*(times/weibullScale)^(weibullShape-1))
  }
  
  theoreticalHazard = sapply(progression_speeds, function(progression_speed){
    getBaselineHazard(weibullScale = weibullScales[progression_speed], 
                      weibullShape = weibullShapes[progression_speed], times = timePoint)
  })
  
  unscaledWeights = sapply(progression_speeds, function(progression_speed){
    weibullScale = weibullScales[progression_speed]
    weibullShape = weibullShapes[progression_speed]
    return(exp(-(timePoint/weibullScale)^weibullShape))
  })
  weights = unscaledWeights/sum(unscaledWeights)
  
  return(sum(theoreticalHazard * weights))
}

getSplineLogBaselineHazard = function(jointModelData, timePoint){
  splineDesign(jointModelData$mvJoint_dre_psa_simDs$control$knots, timePoint, 
               ord = jointModelData$mvJoint_dre_psa_simDs$control$ordSpline, outer.ok = T) %*% jointModelData$mvJoint_dre_psa_simDs$statistics$postMeans$Bs_gammas
}

getSemiParametricLogBaselineHazard = function(jointModelData, semiParamCoeffs, knots, timePoint){
  semiParamHazard = semiParamCoeffs[max(which(knots<=timePoint))]
}

times = seq(0.05, 10, length.out = 500)

savedFiles = list.files(path = "C:/Users/838035/Google Drive/PhD/src/prias/Rdata/decision_analytic/Simulation/Mixed/", full.names = T, pattern = "normal.Rdata")
logHazardMatrix = matrix(ncol=length(times), nrow=length(savedFiles))
jmFitLogHazardMatrix = matrix(ncol=length(times), nrow=length(savedFiles))

for(i in 1:length(savedFiles)){
  rm(jointModelData)
  print(paste("Loading the file number:", i))
  load(savedFiles[i])
  
  semiParamCoeffs = c(semiParamBH(jointModelData$survModel_simDs, length.knots = 20)[[1]]$coefficients[-c(1:length(jointModelData$survModel_simDs$coefficients))], 0)
  knots = semiParamBH(jointModelData$survModel_simDs, length.knots = 20)[[2]]  
  logHazardMatrix[i, ] = sapply(times, getSemiParametricLogBaselineHazard, jointModelData=jointModelData, semiParamCoeffs=semiParamCoeffs, knots=knots)
  #logHazardMatrix[i, ] = sapply(times, getSplineLogBaselineHazard, jointModelData=jointModelData)
  
  print(paste("******** End working on file number: ", i, "*******"))
}

#mixtureTheoreticalLogHazard = log(sapply(times, getTheoreticalHazard, jointModelData=jointModelData))
mixtureTheoreticalLogHazard = splineDesign(mvJoint_dre_psa_dre_value_superlight$control$knots, times, 
                 ord = mvJoint_dre_psa_dre_value_superlight$control$ordSpline, outer.ok = T) %*% mvJoint_dre_psa_dre_value_superlight$statistics$postMeans$Bs_gammas

fittedLogHazardMean = (apply((logHazardMatrix), 2, mean, na.rm=T))
fittedLogHazardLower = (apply((logHazardMatrix), 2, quantile, probs=c(0.025), na.rm=T))
fittedLogHazardUpper = (apply((logHazardMatrix), 2, quantile, probs=c(0.975), na.rm=T))

df = data.frame(logHazard = c(mixtureTheoreticalLogHazard, fittedLogHazardMean),
                "Baseline Hazard"=rep(c("Theoretical log (baseline hazard)", 
                                        "Fitted log(baseline hazard)"), each=500),
                confIntLow=rep(fittedLogHazardLower,2),
                confIntHigh=rep(fittedLogHazardUpper,2),
                time = rep(times,2))

ggplot(data=df[df$time > 0 & df$time<=10,]) + 
  geom_line(aes(x=time, y=logHazard, linetype=Baseline.Hazard)) + 
  geom_ribbon(aes(x=time, ymin=(confIntLow), ymax=(confIntHigh)), fill = "grey", alpha=0.4) + 
  scale_x_continuous(breaks = c(0.1, seq(1, 10, 1))) + 
  scale_linetype_manual(values=c("twodash", "solid")) +
  theme(legend.title=element_blank(), legend.position = "top",
        text = element_text(size=11), axis.text=element_text(size=11),
        legend.text = element_text(size=11)) + 
  xlab("Time (years)") + ylab("log (baseline hazard)")

# ggsave(file="report/pers_schedule/biometrics_submission/images/sim_study/baseline_hazard.eps", 
#        width=8.27*0.75, height=8.27*0.75, device=cairo_ps)
