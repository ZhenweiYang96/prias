load("Rdata/decision_analytic/DRE_PSA/mvJoint_dre_psa_dre_value_light.Rdata")
load("Rdata/decision_analytic/cleandata.Rdata")
source("src/decision_analytic/load_lib.R")
source("src/decision_analytic/simulationStudy/simCommonNoAlpha.R")
source("src/decision_analytic/simulationStudy/schedules.R")

MAX_FAIL_TIME = 15
max_cores = 8
dataSetNums = 1:1

lastSeed = 101
progression_type="Mixed"
for(i in dataSetNums){
  #Save RAM
  rm(jointModelData)  
  print(paste("******** Started working on Data Set: ", i, "*******"))
  
  lastSeed = getNextSeed(lastSeed)
  repeat{
    print(paste("Using seed:", lastSeed))
    jointModelData = try(fitJointModelOnNewData(seed = lastSeed, 
                                                nSub = 2500, nSubTraining = 1500, nSubTest = 250,
                                                censStartTime = 25, censEndTime = 25, 
                                                engine="STAN", progression_type = progression_type),T)
    if(inherits(jointModelData, "try-error")){
      print(jointModelData)
      lastSeed = getNextSeed(lastSeed)
      print("Error: trying again")
    }else{
      break
    }
  }

  saveName = paste0("jointModelData_seed_",jointModelData$seed,"_simNr_",i,"_.Rdata")
  save(jointModelData, file = paste0("Rdata/decision_analytic/Simulation/", progression_type, "/", saveName))
}
