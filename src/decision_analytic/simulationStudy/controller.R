load("Rdata/decision_analytic/DRE_PSA/mvJoint_dre_psa_dre_value_light.Rdata")
load("Rdata/decision_analytic/cleandata.Rdata")
source("src/decision_analytic/load_lib.R")
source("src/decision_analytic/simulationStudy/simCommon.R")
source("src/decision_analytic/simulationStudy/schedules.R")

MAX_FAIL_TIME = 15
max_cores = 8
dataSetNums = 1:10

ANNUAL = "Annual"
MONTH_18 = "18 Months"
BIENNIAL = "Biennial"
PRIAS = "PRIAS"
KAPPApt95 = "Dyn. Risk (5%) GR"
KAPPApt90 = "Dyn. Risk (10%) GR"
KAPPApt85 = "Dyn. Risk (15%) GR"
KAPPApt80 = "Dyn. Risk (20%) GR"
KAPPAF1Score = "Dyn. Risk (F1) GR"
methodNames = c(ANNUAL, MONTH_18, BIENNIAL, PRIAS, KAPPApt95, KAPPApt90, KAPPApt85, KAPPApt80, KAPPAF1Score)

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
                                                nSub = 1500, nSubTraining = 1000, nSubTest = 250,
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
  
  scheduleResults = data.frame(P_ID = jointModelData$testData$testDs.id$P_ID,
                                                       progression_speed = jointModelData$testData$testDs.id$progression_speed,
                                                       progression_time = jointModelData$testData$testDs.id$progression_time,
                                                       Age = jointModelData$testData$testDs.id$Age,
                                                       methodName = rep(methodNames, each=nrow(jointModelData$testData$testDs.id)),
                                                       nb = NA, offset=NA)
  
  #1st we do the annual schedule
  print("Running Annual, 18 month and Biennial schedules")
  scheduleResults[scheduleResults$methodName == ANNUAL, c("nb", "offset")] = runFixedSchedule(jointModelData$testData$testDs.id, biopsyTimes = seq(0, 30, 1))
  scheduleResults[scheduleResults$methodName == MONTH_18, c("nb", "offset")] = runFixedSchedule(jointModelData$testData$testDs.id, biopsyTimes = seq(0, 30, 1.5))
  scheduleResults[scheduleResults$methodName == BIENNIAL, c("nb", "offset")] = runFixedSchedule(jointModelData$testData$testDs.id, biopsyTimes = seq(0, 30, 2))
  print("Done running Annual 18 month and Biennial schedules")
  
  #Then we do the PRIAS schedule
  print("Running PRIAS schedule")
  scheduleResults[scheduleResults$methodName == PRIAS, c("nb", "offset")] = runPRIASSchedule(jointModelData$testData$testDs.id, jointModelData$testData$testDs)
  print("Done running PRIAS schedule")
  
  #Then we do the schedule with Dyn. Risk of GR
  for(riskMethodName in c(KAPPApt95, KAPPApt90, KAPPApt85, KAPPApt80, KAPPAF1Score)){
    print(paste("Running",riskMethodName,"schedule"))
    scheduleResults[scheduleResults$methodName == riskMethodName, c("nb", "offset")] = runDynRiskGRSchedule(jointModelData, riskMethodName)
    print(paste("Done running",riskMethodName,"schedule"))
  }
  
  print(paste("********* Saving the results ******"))
  
  jointModelData$testData$scheduleResults = scheduleResults
  saveName = paste0("jointModelData_seed_",jointModelData$seed,"_simNr_",i,"_.Rdata")
  save(jointModelData, file = paste0("Rdata/decision_analytic/Simulation/", progression_type, "/", saveName))
  rm(scheduleResults)
}
