setwd("C:/Users/anirudhtomer/Google Drive/PhD/src/prias/")
args = commandArgs(trailingOnly=TRUE)

seed = 2039
max_cores = 6

library(JMbayes)
library(splines)
library(survival)
library(MASS)

load("Rdata/gap3/PRIAS_2019/cleandata.Rdata")
source("src/lastpaper/simulation/simCommon.R")

nSub = 1000
psaErrorDist = "t3"
MAX_FAIL_TIME = 10

print(paste("Using seed: ", seed, " and max_cores", max_cores))
set.seed(seed)
newData = generateSimulationData(nSub = nSub, psaErrorDist = psaErrorDist, bNames=bNames)
print('Generated data, now gonna fit the model')
jointModelData = fitJointModelOnNewData(simDs = newData$simDs, simDs.id=newData$simDs.id,
                                        timesPerSubject = newData$timesPerSubject,
                                        nSubTraining = 750, nSubTest = 250, 
                                        censStartTime = 1, censEndTime = 15)
jointModelData$seed = seed

print("********* Saving the model ******")
saveName = paste0("jointModelData_seed_",jointModelData$seed, "_",psaErrorDist, ".Rdata")
save(jointModelData, file = paste0("Rdata/lastpaper/simulation/full/", saveName))

jointModelData$mvglmer_dre_psa_simDs = NULL
jointModelData$mvglmer_fitting_time = NULL
jointModelData$mvjoint_fitting_time = NULL
jointModelData$survModel_simDs = NULL
jointModelData$mvJoint_dre_psa_simDs$mcmc$b = NULL
jointModelData$mvJoint_dre_psa_simDs$call = NULL
jointModelData$mvJoint_dre_psa_simDs$weights = NULL
jointModelData$mvJoint_dre_psa_simDs$mcmc_info = NULL
jointModelData$mvJoint_dre_psa_simDs$model_info = NULL
jointModelData$mvJoint_dre_psa_simDs$control = list(knots = jointModelData$mvJoint_dre_psa_simDs$control$knots,
                                                    ordSpline = jointModelData$mvJoint_dre_psa_simDs$control$ordSpline)

save(jointModelData, file = paste0("Rdata/lastpaper/simulation/light/", saveName))
