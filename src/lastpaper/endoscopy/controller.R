args = commandArgs(trailingOnly=TRUE)
seed = as.numeric(args[1])
max_cores = as.numeric(args[2])

load("Rdata/lastpaper/endoscopy/data.RData")
source("src/lastpaper/endoscopy/simCommon.R")

MAX_FAIL_TIME = 14
n_sub = 1000
n_training = 750

FIXED = "Fixed"
thresholds = 1 - c(1.5,2,2.5,3)/100
min_gaps = c(0.5,1)

print(paste("Using seed:", seed))
set.seed(seed)
new_data = generateBaselineData(n_sub = n_sub)

joint_model_data = try(fitJointModel(new_data$patient_data, 
                                     n_training, cens_start_time = 1, cens_end_time = 30),T)
joint_model_data$seed = last_seed
joint_model_data$b = new_data$b

saveName = paste0("joint_model_data_seed_",joint_model_data$seed,"_simNr_",i, ".Rdata")
save(joint_model_data, file = paste0("Rdata/simulation/", saveName))

new_method_names = c(FIXED, sapply(thresholds, function(x){
  paste0("Risk_",x,"_gap_", min_gaps)
}))

schedule_results = do.call(rbind, replicate(length(new_method_names), 
                                            joint_model_data$test_data[!duplicated(joint_model_data$test_data$id),], 
                                            simplify = F))
schedule_results = schedule_results[order(schedule_results$id, decreasing = F),]

schedule_results$methodName = new_method_names
schedule_results$nb = schedule_results$offset = NA

# Then we do the FIXED schedule
print("Running Fixed schedule")
schedule_results[schedule_results$methodName == FIXED, c("nb", "offset")] = runFixedSchedule(joint_model_data$test_data)
print("Done running Fixed schedule")

#Then we do the schedule with Dyn. Risk of GR
for(min_gap in min_gaps){ 
  #Then we do the schedule with Dyn. Risk of GR
  for(threshold in thresholds){
    riskMethodName = paste0("Risk_",threshold,"_gap_", min_gap)
    
    print(paste("Running",riskMethodName,"schedule"))
    schedule_results[schedule_results$methodName == riskMethodName, c("nb", "offset")] = runRiskBasedSchedule(joint_model_data$fitted_jm,
                                                                                                              joint_model_data$test_data, 
                                                                                                              threshold, min_gap)
    print(paste("Done running",riskMethodName,"schedule"))
  }
}

# print(paste("********* Saving the results ******"))
joint_model_data$schedule_results = schedule_results
save(joint_model_data, file = paste0("Rdata/simulation/", saveName))
rm(schedule_results)
rm(joint_model_data)
