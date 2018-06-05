load("Rdata/decision_analytic/PSA_Only/mvJoint_psa_light.Rdata")
load("Rdata/decision_analytic/cleandata.Rdata")
source("src/decision_analytic/load_lib.R")
source("src/decision_analytic/simulationStudy/PSA_Only/simCommonPSA.R")

MAX_FAIL_TIME = 15
max_cores = 4
dataSetNums = 1:10

methods = c("expectedFailureTime", "medianFailureTime","youden", "f1score") 
#"accuracy", "f1score")

lastSeed = 401
progression_type="Mixed"
for(i in dataSetNums){
  #Save RAM
  rm(jointModelData)  
  print(paste("******** Started working on Data Set: ", i, "*******"))
  
  lastSeed = getNextSeed(lastSeed)
  repeat{
    print(paste("Using seed:", lastSeed))
    jointModelData = try(fitJointModelOnNewData(seed = lastSeed, 
                                                nSub = 1500, nSubTraining = 1000, nSubTest = 250,
                                                censStartTime = 25, censEndTime = 25, 
                                                engine="JAGS", progression_type = progression_type),T)
    if(inherits(jointModelData, "try-error")){
      print(jointModelData)
      lastSeed = getNextSeed(lastSeed)
      print("Error: trying again")
    }else{
      break
    }
  }
  
  print(paste("********* Done working on Data Set"))
  
  saveName = paste0("jointModelData_seed_",jointModelData$seed,"_simNr_",i,"_.Rdata")
  save(jointModelData, file = paste0("Rdata/decision_analytic/Simulation/", progression_type, "/", saveName))
}
