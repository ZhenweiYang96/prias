library(JMbayes)
library(splines)
library(doParallel)

#First load a simulated jointmodel data object
files = list.files("C:/Users/838035/prias/Rdata/mdp/final_res", full.names = T)

delay_list = vector("list", length(files))
for(file_num in 1:length(files)){
  print(paste("File: ", file_num))
  load(files[file_num])
  
  fitted_JM = jointModelData$mvJoint_dre_psa_simDs
  jointModelData$mvJoint_dre_psa_simDs = NULL
  patient_data_list = split(jointModelData$testData$testDs, 
                            f = jointModelData$testData$testDs$P_ID)
  
  ct = makeCluster(4)
  registerDoParallel(ct)
  #First we check error when G > 7 is marked as G < 7
  
  trueG7_markedG6 = foreach(p=1:length(patient_data_list), 
                            .packages = c("splines", "JMbayes")) %dopar%{
                              
                              patient_df = patient_data_list[[p]]
                              prog_time = patient_df$progression_time[1]
                              
                              startIndex = max(which(patient_df$visitTimeYears <= prog_time)) + 1
                              survProbs=c()
                              if(startIndex < nrow(patient_df) && p!=83){
                                for(i in startIndex:nrow(patient_df)){
                                  curTime = patient_df$visitTimeYears[i]
                                  temp = survfitJM(object = fitted_JM, newdata = patient_df[1:i,], idVar="P_ID", 
                                                   last.time=prog_time, survTimes = curTime)
                                  survProbs = c(survProbs, temp$summaries[[1]][1, "Mean"])
                                }
                                names(survProbs) = patient_df$visitTimeYears[startIndex:nrow(patient_df)]
                              }
                              
                              return(survProbs)
                            }
  
  stopCluster(ct)
  
  cutoff = 0.85
  delay_pt15 = sapply(1:length(trueG7_markedG6), function(p){
    patient_df = patient_data_list[[p]]
    prog_time = patient_df$progression_time[1]
    
    survProbs = trueG7_markedG6[[p]]
    if(!is.null(survProbs)){
      survProbs = survProbs - cutoff
      detectTime = min(as.numeric(names(survProbs))[survProbs <= 0], 10, na.rm = T)
      return(detectTime - prog_time)
    }
    return(NA)
  })
  
  cutoff = 0.9
  delay_pt1 = sapply(1:length(trueG7_markedG6), function(p){
    patient_df = patient_data_list[[p]]
    prog_time = patient_df$progression_time[1]
    
    survProbs = trueG7_markedG6[[p]]
    if(!is.null(survProbs)){
      survProbs = survProbs - cutoff
      detectTime = min(as.numeric(names(survProbs))[survProbs <= 0], 10, na.rm = T)
      return(detectTime - prog_time)
    }
    return(NA)
  })
  
  cutoff = 0.95
  delay_pt05 = sapply(1:length(trueG7_markedG6), function(p){
    patient_df = patient_data_list[[p]]
    prog_time = patient_df$progression_time[1]
    
    survProbs = trueG7_markedG6[[p]]
    if(!is.null(survProbs)){
      survProbs = survProbs - cutoff
      detectTime = min(as.numeric(names(survProbs))[survProbs <= 0], 10, na.rm = T)
      return(detectTime - prog_time)
    }
    return(NA)
  })
  
  delay_list[[file_num]] = data.frame(prog_time=rep(jointModelData$testData$testDs.id$progression_time,3),
                                      delay=c(delay_pt15, delay_pt1, delay_pt05),
                                      cutoff=rep(c(0.15, 0.1, 0.05), each=length(delay_pt05)))
}

delay_df = do.call('rbind', delay_list)
temp = delay_df[!is.na(delay_df$delay) & delay_df$cutoff==0.1,]
by(INDICES = temp$prog_time<3.5, temp$delay, summary)
ggplot(data=temp) + geom_boxplot(aes(x=factor(cutoff), y=delay), outlier.shape = NA) + 
  facet_grid(.~(prog_time>3.5))

