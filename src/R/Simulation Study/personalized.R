#The weibull scale and shape matters
minVisits = 5
patientDsList = split(testDs, testDs$P_ID)
# 
# res = foreach(i=1:nrow(testDs.id),
#         .packages = c("splines", "JMbayes")) %dopar%{
          
for(i in 1:50){
  patientDs_i = patientDsList[[i]]
  
  patientDs_i$expectedFailureTime = NA
  patientDs_i$survTime90 = NA
  patientDs_i$survTime80 = NA
  
  biopsy_times = list(expectedFailureTime = c(), survTime90 = c(), survTime80 = c())
  
  for(j in minVisits:(timesPerSubject-1)){
    persTestDs = patientDs_i[patientDs_i$visitNumber <= j,]
    
    upperLimitIntegral = 100
    repeat{
      condFailTime = try(expectedCondFailureTime(persTestDs, upperLimitIntegral = upperLimitIntegral))
      
      if(!inherits(condFailTime, "try-error")){
        patientDs_i$expectedFailureTime[j] = condFailTime
        break
      }else{
        upperLimitIntegral = upperLimitIntegral - 10
      }
    }
    patientDs_i$survTime90[j] = pDynSurvTime(survProb = 0.9, persTestDs)
    patientDs_i$survTime80[j] = pDynSurvTime(survProb = 0.8, persTestDs)
    
    if(patientDs_i$expectedFailureTime[j] < patientDs_i$visitTimeYears[j + 1]){
      biopsy_times$expectedFailureTime = c(biopsy_times$expectedFailureTime, patientDs_i$expectedFailureTime[j])
    }
    
    if(patientDs_i$survTime90[j] < patientDs_i$visitTimeYears[j + 1]){
      biopsy_times$survTime90 = c(biopsy_times$survTime90, patientDs_i$survTime90[j])
    }
    
    if(patientDs_i$survTime80[j] < patientDs_i$visitTimeYears[j + 1]){
      biopsy_times$survTime80 = c(biopsy_times$survTime80, patientDs_i$survTime80[j])
    }
  }
  
  res[[i]]=list(patientDs_i, biopsy_times)
}

save.image("Rdata/Gleason as event/simStudy1000_scale.Rdata")

#Store original results as they are
orig_res = res

###########################################################
#Approach #1 do biopsy as given by the estimates from the above model
############################################################
approach1res = foreach(i=1:length(res), .combine = "rbind")%do%{
  trueProgressionTime = testDs.id[i,]$progression_time
  
  biopsy_times = res[[i]][[2]]
  #Expected value
  nb1 = which(biopsy_times$expectedFailureTime >= trueProgressionTime)[1]
  nb2 = which(biopsy_times$survTime90 >= trueProgressionTime)[1]
  nb3 = which(biopsy_times$survTime80 >= trueProgressionTime)[1]
  
  biopsytimeOffset1 = biopsytimeOffset2 = biopsytimeOffset3 = NA
  
  if(!is.na(nb1)){
    biopsytimeOffset1 = biopsy_times$expectedFailureTime[nb1] - trueProgressionTime
  }
  
  if(!is.na(nb2)){
    biopsytimeOffset2 = biopsy_times$survTime90[nb2] - trueProgressionTime
  }
  
  if(!is.na(nb3)){
    biopsytimeOffset3 = biopsy_times$survTime80[nb2] - trueProgressionTime
  }
  
  c(nb1, nb2, nb3, biopsytimeOffset1, biopsytimeOffset2, biopsytimeOffset3)
}

#############################################################
#Approach #2 do biopsy if the failure time is less than an year then do biopsy at that proposed time
############################################################
approach2res = foreach(i=1:length(res), .combine = "rbind")%do%{
  trueProgressionTime = testDs.id[i,]$progression_time
  
  biopsy_times = res[[i]][[2]]
  patientDs_i = res[[i]][[1]]
  
  nb1 = nb2 = nb3 = nb4= 0
  biopsytimeOffset1 = biopsytimeOffset2 = biopsytimeOffset3 = biopsytimeOffset4 = NA
  
  #Expected value
  last_biopsy_time = -Inf
  for(j in minVisits:(timesPerSubject-1)){
    biopsy_gap_needed = (patientDs_i$visitTimeYears[j] - last_biopsy_time) < 1
      
    if(biopsy_gap_needed==FALSE & (patientDs_i$expectedFailureTime[j] - patientDs_i$visitTimeYears[j]) <= 1){
      nb1 = nb1 + 1
      biopsytimeOffset1 = patientDs_i$expectedFailureTime[j] - trueProgressionTime
      last_biopsy_time = patientDs_i$expectedFailureTime[j]
        
      if(biopsytimeOffset1 >= 0){
        break
      }
    }
  }
  
  #survtime90 value
  last_biopsy_time = -Inf
  for(j in minVisits:(timesPerSubject-1)){
    biopsy_gap_needed = (patientDs_i$visitTimeYears[j] - last_biopsy_time) < 1
    
    if(biopsy_gap_needed==FALSE & (patientDs_i$survTime90[j] - patientDs_i$visitTimeYears[j]) <= 1){
      nb2 = nb2 + 1
      biopsytimeOffset2 = patientDs_i$survTime90[j] - trueProgressionTime
      last_biopsy_time = patientDs_i$survTime90[j]
      
      if(biopsytimeOffset2 >= 0){
        break
      }
    }
  }
  
  #survtime80 value
  last_biopsy_time = -Inf
  for(j in minVisits:(timesPerSubject-1)){
    biopsy_gap_needed = (patientDs_i$visitTimeYears[j] - last_biopsy_time) < 1
    
    if(biopsy_gap_needed==FALSE & (patientDs_i$survTime80[j] - patientDs_i$visitTimeYears[j]) <= 1){
      nb3 = nb3 + 1
      biopsytimeOffset3 = patientDs_i$survTime80[j] - trueProgressionTime
      last_biopsy_time = patientDs_i$survTime80[j]
      
      if(biopsytimeOffset3 >= 0){
        break
      }
    }
  }
  
  #Fixed schedule
  if(trueProgressionTime<=1){
    nb4 = 1
    biopsytimeOffset4 = 1 - trueProgressionTime
  }else if(trueProgressionTime<=4){
    nb4 = 2
    biopsytimeOffset4 = 4 - trueProgressionTime
  }else if(trueProgressionTime<=7){
    nb4 = 3
    biopsytimeOffset4 = 7 - trueProgressionTime
  }else if(trueProgressionTime<=10){
    nb4 = 4
    biopsytimeOffset4 = 10 - trueProgressionTime
  }else if(trueProgressionTime<=15){
    nb4 = 5
    biopsytimeOffset4 = 15 - trueProgressionTime
  }else if(trueProgressionTime<=20){
    nb4 = 6
    biopsytimeOffset4 = 20 - trueProgressionTime
  }
  
  c(nb1, nb2, nb3, nb4, biopsytimeOffset1, biopsytimeOffset2, biopsytimeOffset3, biopsytimeOffset4)
}

approach2res = approach2res[testDs.id$progression_time>1 & !is.na(approach2res[,5]) & !is.na(approach2res[,6]) & !is.na(approach2res[,7]) & approach2res[,5]>0 & approach2res[,6]>0 & approach2res[,7]>0,]
apply(approach2res, MARGIN = 2, mean, na.rm=T)
