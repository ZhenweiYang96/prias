getTheoreticalHazard = function(timePoint, progression_speeds){
  
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