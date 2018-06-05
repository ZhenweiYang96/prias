library(MASS)
library(splines)

seeds = 6000:6020
Dalphas1 = rep(NA, length(seeds))
Dalphas2 = rep(NA, length(seeds))
nSub= 500
load("Rdata/decision_analytic/testingDimSim_quadeff.RData")

# parameters for the linear mixed effects model
betas <- jmFit_both$coefficients$betas
sigma.y <- jmFit_both$coefficients$sigma

# parameters for the survival model
gammas <- jmFit_both$coefficients$gammas[-1]
alphas <- c(jmFit_both$coefficients$alpha, jmFit_both$coefficients$Dalpha)

D <- jmFit_both$coefficients$D

for(seedCounter in 1:length(seeds)){
  
  print(seedCounter)
  set.seed(seeds[seedCounter])

  visitTimes = seq(0, 10, 0.05)
  id = rep(1:nSub, each=length(visitTimes))
  times <- rep(visitTimes, nSub)
  
  drug = rep(sample(c("D-penicil","placebo"), replace = T, size = nSub), each=length(visitTimes))
  age = rep(rnorm(nSub, mean=mean(pbc2.id$age), sd = sd(pbc2.id$age)), each=length(visitTimes))
  
  DF <- data.frame(id=id, year = times, drug = factor(drug, labels = c("placebo", "D-penicil")), age=age)
  DF.id = DF[!duplicated(DF$id),]
  
  fixedFormula = ~ 1 + drug + age + I(year-3) + I((year-3)^2)
  randomFormula = ~ 1 + I(year-3) + I((year-3)^2)
  
  fixedSlopeFormula = ~ 1 + I(2*(year-3))
  randomSlopeFormula = ~ 1 + I(2*(year-3))
  
  X <- model.matrix(fixedFormula, data = DF)
  Z <- model.matrix(randomFormula, data = DF)
  
  # design matrix for the survival model
  W <- model.matrix(~1+drug+age, DF.id)[,-1]
  
  ################################################
  
  #simulate random effects
  b <- mvrnorm(nSub, rep(0, nrow(D)), D)
  # simulate longitudinal responses
  eta.y <- as.vector(X %*% betas + rowSums(Z * b[id, ])) # linear predictor
  y <- rnorm(nSub * length(visitTimes), eta.y, sigma.y)
  DF$y = y
  # simulate event times
  # weibullShape = 6
  # weibullScale = 11
  
  #weibullShapes = c(5,6)
  #weibullScales = c(8,11)
  
  weibullShapes = c(4,6)
  weibullScales = c(7,11)
  
  eta.t <- as.vector(W %*% gammas)
  invS <- function (t, u, i) {
    index = sample(1:2,1)
    index = 1
    weibullShape = weibullShapes[index]
    weibullScale = weibullScales[index]
    
    h <- function (s) {
      tempDf = data.frame(year=s, drug=DF.id$drug[i], age=DF.id$age[i])
      XX <- model.matrix(fixedFormula, tempDf)
      ZZ <- model.matrix(randomFormula, tempDf)
      
      XX_slope <- model.matrix(fixedSlopeFormula, tempDf)
      ZZ_slope <- model.matrix(randomSlopeFormula, tempDf)
      
      f1 <- as.vector(XX %*% betas + ZZ %*% b[i, ])
      f2 = as.vector(XX_slope %*% betas[-c(1:3)] + ZZ_slope %*% b[i, -1])
      baselinehazard = ((weibullShape/weibullScale)*(s/weibullScale)^(weibullShape-1))
      baselinehazard * exp(eta.t[i] + cbind(f1,f2) %*% alphas)
    }
    integrate(h, lower = 0, upper = t)$value + log(u)
  }
  u <- runif(nSub)
  DF.id$Time <- numeric(nSub)
  DF.id$event = 1
  for (i in 1:nSub) {
    Up <- 50
    tries <- 5
    Root <- try(uniroot(invS, interval = c(1e-05, Up), u = u[i], i = i)$root, TRUE)
    while(inherits(Root, "try-error") && tries > 0) {
      tries <- tries - 1
      Up <- Up + 200
      Root <- try(uniroot(invS, interval = c(1e-05, Up), u = u[i], i = i)$root, TRUE)
    }
    DF.id$Time[i] <- if (!inherits(Root, "try-error")) Root else NA
  }
  
  summary(DF.id$Time)
  
  na.id <- which(is.na(DF.id$Time))
  DF = DF[!DF$id %in% na.id,]
  DF.id = DF.id[!DF.id$id %in% na.id,]
  DF$Time = rep(DF.id$Time, each=length(visitTimes))
  
  DF = DF[DF$year <= DF$Time, ]
  
  # ################################################
  # 
  # Fit the corresponding joint model using package JM
  mvglmer_training <- mvglmer(list(y~ drug + age + I(year-3) + I((year-3)^2) +
                                     ((I(year-3) + I((year-3)^2)) | id)), data = DF,
                              families = list(gaussian))

  # [2] Fit a Cox model, specifying the baseline covariates to be included in the joint
  # model; you need to set argument 'model' to TRUE.
  CoxFit_training <- coxph(Surv(Time, event) ~ drug + age, data = DF.id, model = TRUE,x = T)

  # [3] The basic joint model is fitted using a call to mvJointModelBayes(), which is very
  # similar to JointModelBayes(), i.e.,
  mvJM_training <- mvJointModelBayes(mvglmer_training, CoxFit_training, timeVar = "year",
                                     Formulas = list("y" = "value",
                                                     "y" = list(fixed = ~ 1 + I(2*(year-3)), random=~1 + I(2*(year-3)),
                                                                indFixed = 4:5, indRandom=2:3, name = "slope")),
                                     control=list(n_cores=8))
  # ##############
  # 
  # lmeFit <- lme(log(serBilir) ~ drug + age + I(year-3) + I((year-3)^2),
  #                random = ~(I(year-3) + I((year-3)^2))|id,
  #                data = pbc2, control = lmeControl(opt = "optim", optimMethod = "L-BFGS-B"))
  # 
  # CoxFit <- coxph(Surv(Time, event) ~ drug + age, data = pbc2.id, model = TRUE,x = T)
  # jmFit_value <- jointModel(lmeFit, CoxFit, timeVar = "year", method = "weibull-PH-aGH")
  # jmFit_both = jointModel(lmeFit, CoxFit, timeVar = "year", method = "weibull-PH-aGH",
  #                          parameterization = "both",
  #                          derivForm = list(fixed = ~ 1 + I(2*(year-3)), random=~1 + I(2*(year-3)),
  #                   indFixed = 4:5, indRandom=2:3))
  
  
  # lmeFit_training <- lme(y~ drug + age + I(year-3) + I((year-3)^2),
  #                        random = ~(I(year-3) + I((year-3)^2))|id,
  #                        data = DF, control = lmeControl(opt = "optim", optimMethod = "L-BFGS-B"))
  # 
  # CoxFit_training <- coxph(Surv(Time, event) ~ drug + age, data = DF.id, model = TRUE,x = T)
  # #jmFit_value_training <- jointModel(lmeFit, CoxFit, timeVar = "year", method = "weibull-PH-aGH")
  # jmFit_both_training = jointModel(lmeFit_training, CoxFit_training, timeVar = "year", method = "weibull-PH-aGH",
  #                                  parameterization = "both",
  #                                  derivForm = list(fixed = ~ 1 + I(2*(year-3)), random=~1 + I(2*(year-3)),
  #                                                   indFixed = 4:5, indRandom=2:3))
  
  Dalphas1[seedCounter] = mvJM_training$statistics$postwMeans$alphas[2]
  Dalphas2[seedCounter] = mvJM_training$statistics$postMeans$alphas[2]
}
