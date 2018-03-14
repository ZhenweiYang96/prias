library(MASS)
library(splines)

nSub <- 1000# number of subjects

source("src/R/common.R")
source("src/R/ascertainment_bias/simCommon.R")

cores = detectCores()

#3 types of weibull scales and shapes
weibullScales = c(4,5,6)
weibullShapes = c(1.5,3,4.5)

#weibullScales = c(5,5,5)
#weibullShapes = c(3,3,3)

simulatedDsList = vector("list", 1)
set.seed(176)
i=1
simulatedDsList[[i]] = generateSimLongitudinalData(nSub)
indices = sample(1:length(weibullScales), size = nSub, replace = T, prob = rep(1/length(weibullScales),length(weibullScales)))
simulatedDsList[[i]]$weibullShape = weibullShapes[indices]
simulatedDsList[[i]]$weibullScale = weibullScales[indices]
    
simulatedDsList[[i]] = try(generateSimJointData(i, nSub, psa_dt_yes = F, intervalCensoring = F))
print(paste("Nr. of training:", nrow(simulatedDsList[[i]]$trainingDs.id), 
              "; Nr. of test:", nrow(simulatedDsList[[i]]$testDs.id)))
  
print("Beginning to fit joint model")
simulatedDsList[[i]]$models = fitJointModelSimDs(simulatedDsList[[i]]$trainingDs.id, simulatedDsList[[i]]$trainingDs, intervalCensoring = F)
print("Joint model fitted")
  
