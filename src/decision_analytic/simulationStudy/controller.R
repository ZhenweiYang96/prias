load("Rdata/decision_analytic/DRE_PSA/mvJoint_dre_psa_dre_value_superlight.Rdata")
source("src/decision_analytic/load_lib.R")
source("src/decision_analytic/simulationStudy/simCommon.R")
source("src/decision_analytic/simulationStudy/schedules.R")

## Weibull shape for constant has to be 1, scale can be any value > 0
# The weibull location tells the minimum failure time.
# weibullScales = c("Fast"=3, "Medium"=5, "Slow"=8)
# weibullShapes = c("Fast"=5, "Medium"=8, "Slow"=14)
# weibullLocations = c("Fast"=0, "Medium"=0, "Slow"=0, "Constant"=0)

psaErrorDist = "normal"
#fail_time_search_upper_limit = 10
#MAX_FAIL_TIME = max(seq(0,fail_time_search_upper_limit,0.1)[which(!is.nan(sapply(seq(0,fail_time_search_upper_limit,0.1), getTheoreticalHazard, progression_speeds=progression_speeds)))])
MAX_FAIL_TIME = 10
MIN_FIXED_ROWS = 2
max_cores = min(8, detectCores())
dataSetNums = 1:10

bNames = c("b_Int_DRE", "b_Slope_DRE", "b_Int_PSA", "b_Slope1_PSA", "b_Slope2_PSA", "b_Slope3_PSA", "b_Slope4_PSA")

ANNUAL = "Annual"
MONTH_18 = "18 Months"
BIENNIAL = "Biennial"
PRIAS = "PRIAS"
EXP_FAIL_TIME = "Mean Time"
MEDIAN_FAIL_TIME = "Median Time"
HYBRID = "Hybrid"

KAPPApt95 = "Risk (5%)"
KAPPApt90 = "Risk (10%)"
KAPPApt85 = "Risk (15%)"
KAPPApt80 = "Risk (20%)"

DYN_RISK_GR = c()
DYN_RISK_GR[KAPPApt95] = 0.95
DYN_RISK_GR[KAPPApt90] = 0.90
DYN_RISK_GR[KAPPApt85] = 0.85
DYN_RISK_GR[KAPPApt80] = 0.80

KAPPAF1Score = "Risk (F1)"
KAPPAYouden = "Risk (J1)"

methodNames = c(ANNUAL, MONTH_18, BIENNIAL, PRIAS, names(DYN_RISK_GR), 
                KAPPAF1Score, EXP_FAIL_TIME, MEDIAN_FAIL_TIME, HYBRID)
                #KAPPAYouden)
                #"Try1", "Try2", "Try3", "Try4", "Try5", "Try6","Try7", "Try8", "Try9", "Try10")

nSub = 1500
lastSeed = 200
for(i in dataSetNums){
  #Save RAM
  rm(jointModelData)  
  print(paste("******** Started working on Data Set: ", i, "*******"))
  
  lastSeed = getNextSeed(lastSeed)
  repeat{
    print(paste("Using seed:", lastSeed))
    newData = try(generateSimulationData(seed=lastSeed, nSub = nSub, psaErrorDist = psaErrorDist, bNames=bNames),T)
    if(inherits(newData, "try-error")){
      print(newData)
      lastSeed = getNextSeed(lastSeed)
      print("Error generating simulation data: trying again")
    }else{
      jointModelData = try(fitJointModelOnNewData(seed = lastSeed, simDs = newData$simDs, simDs.id=newData$simDs.id,
                                                  timesPerSubject = newData$timesPerSubject,
                                              nSubTraining = 1000, nSubTest = 250, mvglmer_iter = 1000, 
                                              censStartTime = 30, censEndTime = 30,engine = "STAN"),T)
      if(inherits(jointModelData, "try-error")){
        print(jointModelData)
        lastSeed = getNextSeed(lastSeed)
        print("Error fitting JM: trying again")
      }else{
        break
      }
    }
  }
  
  scheduleResults = data.frame(P_ID = jointModelData$testData$testDs.id$P_ID,
                                                       progression_time = jointModelData$testData$testDs.id$progression_time,
                                                       Age = jointModelData$testData$testDs.id$Age,
                                                       methodName = rep(methodNames, each=nrow(jointModelData$testData$testDs.id)),
                                                       nb = NA, offset=NA)

  #1st we do the annual schedule
  print("Running Annual, 18 month and Biennial schedules")
  scheduleResults[scheduleResults$methodName == ANNUAL, c("nb", "offset")] = runFixedSchedule(jointModelData$testData$testDs.id, biopsyTimes = seq(1, 10, 1))
  scheduleResults[scheduleResults$methodName == MONTH_18, c("nb", "offset")] = runFixedSchedule(jointModelData$testData$testDs.id, biopsyTimes = seq(1.5, 10, 1.5))
  scheduleResults[scheduleResults$methodName == BIENNIAL, c("nb", "offset")] = runFixedSchedule(jointModelData$testData$testDs.id, biopsyTimes = seq(2, 10, 2))
  print("Done running Annual 18 month and Biennial schedules")

  #Then we do the PRIAS schedule
  print("Running PRIAS schedule")
  scheduleResults[scheduleResults$methodName == PRIAS, c("nb", "offset")] = runPRIASSchedule(jointModelData$testData$testDs.id, jointModelData$testData$testDs)
  print("Done running PRIAS schedule")

  #Then we do the schedule with median failure time
  print("Running median failure time schedule")
  scheduleResults[scheduleResults$methodName == MEDIAN_FAIL_TIME, c("nb", "offset")] = runMedianFailureTimeSchedule(jointModelData)
  print("Done running median failure time schedule")

  #Then we do the schedule with expected failure time
  print("Running expected failure time schedule")
  scheduleResults[scheduleResults$methodName == EXP_FAIL_TIME, c("nb", "offset")] = runExpectedFailureTimeSchedule(jointModelData)
  print("Done running expected failure time schedule")
  
  #Then we do the schedule with Dyn. Risk of GR
  for(riskMethodName in c(KAPPApt95, KAPPApt90, KAPPApt85, KAPPApt80)){
    print(paste("Running",riskMethodName,"schedule"))
    scheduleResults[scheduleResults$methodName == riskMethodName, c("nb", "offset")] = runFixedRiskGRSchedule(jointModelData, DYN_RISK_GR[riskMethodName])
    print(paste("Done running",riskMethodName,"schedule"))
  }
  
  print("Running F1 score schedule")
  scheduleResults[scheduleResults$methodName == KAPPAF1Score, c("nb", "offset")] = runDynRiskGRSchedule(jointModelData, "F1score")
  print("Done running F1 score schedule")
  
  print(paste("********* Saving the results ******"))

  jointModelData$testData$scheduleResults = scheduleResults
  saveName = paste0("jointModelData_seed_",jointModelData$seed,"_simNr_",i,"_",psaErrorDist, ".Rdata")
  save(jointModelData, file = paste0("Rdata/decision_analytic/Simulation/", saveName))
  rm(scheduleResults)
}
