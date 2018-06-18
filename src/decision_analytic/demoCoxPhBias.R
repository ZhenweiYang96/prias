library(MASS)
library(ggplot2)
library(ggpubr)

########################################
# first I will define some functions needed for comparing theoretical and fitted hazard
########################################
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


######################
set.seed(2018)
##########################################
# Step 1: Decide the progression speed for the population using the variable progression_speeds
# Progression speed can be "Fast", "Medium", "Slow", "Constant"
# Weibull shape for constant has to be 1, scale can be any value > 0
##########################################
progression_speeds = c("Fast", "Medium", "Slow")
#progression_speeds = c("Fast")

weibullScales = c("Fast"=3, "Medium"=5, "Slow"=8, "Constant"=1)
weibullShapes = c("Fast"=5, "Medium"=8, "Slow"=14, "Constant"=1)

#############################################################################
# The weibull parameterization I use lead to the following formula for hazard
# https://en.wikipedia.org/wiki/Weibull_distribution#Cumulative_distribution_function
# getBaselineHazard = function(weibullScale, weibullShape, times){
#   return((weibullShape/weibullScale)*(times/weibullScale)^(weibullShape-1))
# }
#############################################################################

###########################################################
# Step 2: Generate patients using PRIAS data
###########################################################
nSub = 1000

#my data set has variables, Age and log_odds_high_dre at year 1.
dataSet.id = data.frame("P_ID"=1:nSub)
dataSet.id$progression_speed = sample(progression_speeds, nSub, replace = T)

dataSet.id$Age = rnorm(n = nSub, mean=70, sd=7)

#######################################################
#Now to add a variable (not time dependent): log_odds_high_dre at year 1
#I am using the parameters from the fitted PRIAS object
#######################################################
dataSet.id$visitTimeYears = 1
b <- mvrnorm(nSub, mu = c(0,0), Sigma = matrix(c(7.292, -0.458, -0.458, 1.199),ncol = 2, nrow=2))
fixedDREFormula = ~ 1 + I(Age - 70) +  I((Age - 70)^2) + visitTimeYears
randomDREFormula = ~ 1 + visitTimeYears
betas_dre = c(-4.117, 0.0512000128, -0.0007820254, -0.5251386778)

dataSet.id$log_odds_high_dre = model.matrix(fixedDREFormula, data=dataSet.id) %*% betas_dre + rowSums(model.matrix(randomDREFormula, data=dataSet.id) * b)
dataSet.id = dataSet.id[, c("P_ID", "Age", "log_odds_high_dre", "progression_speed")]

###################################################
# calculate wGamma part of the PH relative risk model
# covariates in survival model are: age, quadratic age, and log_odds_high_dre at year 1 as calculated above
####################################################
survivalFormula = ~ 0 + I(Age - 70) + I((Age - 70)^2) + log_odds_high_dre
gammas = c(0.01055815, -0.00132913, 0.15216027)
dataSet.id$wGamma = model.matrix(survivalFormula, data = dataSet.id) %*% gammas

########################################################
# Generate event times
########################################################
# dataSet.id$progression_time = sapply(dataSet.id$P_ID, FUN = function(patientId){
#   survProb = runif(n = 1, 0, 1)
# 
#   wGamma = dataSet.id$wGamma[patientId]
#   progression_speed = dataSet.id$progression_speed[patientId]
# 
#   weibullScale = weibullScales[progression_speed]
#   weibullShape = weibullShapes[progression_speed]
# 
#   weibullScale * ((-log(survProb)/exp(wGamma))^(1/weibullShape))
# })

dataSet.id$progression_time = sapply(dataSet.id$P_ID, FUN = function(patientId){
  invSurvival <- function (t, u, patientId) {
    log(u) + integrate(hazardFunc, lower = 0, upper = t, patientId)$value
  }

  hazardFunc = function (s, patientId) {
    wGamma = dataSet.id$wGamma[patientId]

    baselinehazard_s = sapply(s, getTheoreticalHazard, progression_speeds=progression_speeds)
    baselinehazard_s * exp(wGamma)
  }

  survProb = runif(n = 1, 0, 1)
  Low = 1e-05
  Up <- 11
  Root <- try(uniroot(invSurvival, interval = c(Low, Up),
                        u = survProb, patientId)$root, TRUE)
  if(inherits(Root, "try-error")){
    print(Root)
    return(NA)
  }else{
    return(Root)
  }
})

dataSet.id$progressed = 1

########################################
# Fit a CoxPH model
########################################
relative_risk_fitted = coxph(Surv(progression_time, progressed) ~ I(Age - 70) +  I((Age - 70)^2) +  log_odds_high_dre,
      data = dataSet.id, x=T, model=T)

summary(relative_risk_fitted)

#Time points at which to show fitted hazard
times = seq(min(dataSet.id$progression_time), max(dataSet.id$progression_time), length.out = 500)

#Theoretical Baseline hazard
theoreticalLogMixtureBaselineHazard = log(sapply(times, getTheoreticalHazard, progression_speeds=progression_speeds))

componentLogBaselineHazards = sapply(progression_speeds, FUN = function(progression_speed){
  log(sapply(times, getTheoreticalHazard, progression_speeds=progression_speed))
})

colnames(componentLogBaselineHazards) = progression_speeds

semiParamCoeffs = semiParamBH(relative_risk_fitted, length.knots = 20)[[1]]$coefficients[-c(1:length(relative_risk_fitted$coefficients))]
knots = semiParamBH(relative_risk_fitted, length.knots = 20)[[2]]
fittedLogBaselineHazard = sapply(times, FUN = function(timePoint, semiParamCoeffs, knots){
  semiParamHazard = semiParamCoeffs[max(which(knots<=timePoint))]
}, semiParamCoeffs=semiParamCoeffs, knots=knots)

plotDf = data.frame(logHazard = c(theoreticalLogMixtureBaselineHazard, fittedLogBaselineHazard),
                "Baseline Hazard"=rep(c("Theoretical log (mixture baseline hazard)", 
                                        "Fitted log(baseline hazard)"), each=length(times)),
                time = rep(times,2))
for(progression_speed in progression_speeds){
  plotDf = rbind(plotDf, data.frame(logHazard = componentLogBaselineHazards[, progression_speed],
                                    "Baseline Hazard"=rep(paste0("Theoretical log (", progression_speed, " baseline hazard)"), length(times)),
                                    time = times))
}

p1 = ggplot(data=plotDf[1:1000,]) + 
  geom_line(aes(x=time, y=logHazard, linetype=Baseline.Hazard, color=Baseline.Hazard)) + 
  scale_linetype_manual(values=c("twodash", "solid", "dotted", "dotted", "dotted")) +
  theme(legend.title=element_blank(), legend.position = "top",
        text = element_text(size=11), axis.text=element_text(size=11),
        legend.text = element_text(size=11)) + 
  xlab("Progression Time (years)") + ylab("log (baseline hazard)")
print(p1)

p2 = ggplot(data=dataSet.id) + geom_density(aes(x=progression_time, fill=progression_speed, alpha=0.1)) +
  theme(legend.title=element_blank(), legend.position = "bottom",
        text = element_text(size=11), axis.text=element_text(size=11),
        legend.text = element_text(size=11)) + 
  xlab("Progression Time (years)")

ggarrange(p1, p2, ncol = 1, nrow = 2)
