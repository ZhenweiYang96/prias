load("Rdata/decision_analytic/PSA_Only/mvJoint_psa_superlight.Rdata")
source("src/decision_analytic/load_lib.R")
source("src/decision_analytic/simulationStudy/PSA_Only/simCommonPSA.R")

## Weibull shape for constant has to be 1, scale can be any value > 0
# The weibull location tells the minimum failure time.
# weibullScales = c("Fast"=3, "Medium"=5, "Slow"=8)
# weibullShapes = c("Fast"=5, "Medium"=8, "Slow"=14)
# weibullLocations = c("Fast"=0, "Medium"=0, "Slow"=0, "Constant"=0)

weibullScales = c("Fast"=2.25, "Medium"=8.25, "Slow"=10, "Constant"=20)
weibullShapes = c("Fast"=2, "Medium"=3.0, "Slow"=14, "Constant"=1)
weibullLocations = c("Fast"=0, "Medium"=0, "Slow"=0, "Constant"=0)

progression_speeds=c("Fast", "Medium", "Slow")
#progression_speeds = c("Fast", "Medium", "Medium")
subFolderName = "Mixed"

psaErrorDist = "normal"
fail_time_search_upper_limit = 10
MAX_FAIL_TIME = max(seq(0,fail_time_search_upper_limit,0.1)[which(!is.nan(sapply(seq(0,fail_time_search_upper_limit,0.1), getTheoreticalHazard, progression_speeds=progression_speeds)))])
MIN_FIXED_ROWS = 2
max_cores = min(8, detectCores())
dataSetNums = 1:10

bNames = c("b_Int_PSA", "b_Slope1_PSA", "b_Slope2_PSA", "b_Slope3_PSA", "b_Slope4_PSA")

nSub = 1000
lastSeed = 200
for(i in dataSetNums){
  #Save RAM
  rm(jointModelData)  
  print(paste("******** Started working on Data Set: ", i, "*******"))
  
  lastSeed = getNextSeed(lastSeed)
  repeat{
    print(paste("Using seed:", lastSeed))
    newData = try(generateSimulationData(seed=lastSeed, nSub = nSub, psaErrorDist = psaErrorDist, 
                                         progression_speeds = progression_speeds, bNames=bNames),T)
    if(inherits(newData, "try-error")){
      print(newData)
      lastSeed = getNextSeed(lastSeed)
      print("Error generating simulation data: trying again")
    }else{
      jointModelData = try(fitJointModelOnNewData(seed = lastSeed, simDs = newData$simDs, simDs.id=newData$simDs.id,
                                                  timesPerSubject = newData$timesPerSubject,
                                                  nSubTraining = 750, nSubTest = 250, mvglmer_iter = 1000, 
                                                  censStartTime = 30, censEndTime = 30,engine = "JAGS"),T)
      if(inherits(jointModelData, "try-error")){
        print(jointModelData)
        lastSeed = getNextSeed(lastSeed)
        print("Error fitting JM: trying again")
      }else{
        break
      }
    }
  }

  saveName = paste0("jointModelData_seed_",jointModelData$seed,"_simNr_",i,"_",psaErrorDist, ".Rdata")
  save(jointModelData, file = paste0("Rdata/decision_analytic/Simulation/", subFolderName, "/", saveName))
  rm(scheduleResults)
}
