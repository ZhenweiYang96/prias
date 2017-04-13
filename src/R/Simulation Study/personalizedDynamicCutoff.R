dynamicCutOffTimes = seq(generateLongtiudinalTimeBySchedule()[minVisits], 10, 0.1)

#Calculate the prevalence of the disease
kmfit = survfit(Surv(progression_time, progressed)~1, conf.type="log-log", data=prias.id)
kmprob <- stepfun(kmfit$time, c(1, kmfit$surv))
prevalence = kmprob(dynamicCutOffTimes) - kmprob(dynamicCutOffTimes + 1)

rocList = foreach(tstart = dynamicCutOffTimes, .packages = c("splines", "JMbayes")) %dopar%{
  res <- tryCatch({
    rocJM(simJointModel_replaced, trainingDs, Tstart=tstart, Dt=1, idVar = "P_ID")
  }, error=function(e) NULL)
}

youdenDynamicCutoffValues = sapply(1:length(rocList), function(i){
  
  k=0
  repeat{
    rocRes = rocList[[i-k]]
    
    if(is.null(rocRes) | any(is.nan(rocRes$TP)) | any(is.nan(rocRes$FP))){
      k = k + 1
    }else{
      break
    }
  }
  rocRes$thrs[which.max(rocRes$TP - rocRes$FP)]
})
markedNess = sapply(1:length(rocList), function(i){
  
  k=0
  repeat{
    rocRes = rocList[[i-k]]
    
    if(is.null(rocRes) | any(is.nan(rocRes$TP)) | any(is.nan(rocRes$FP))){
      k = k + 1
    }else{
      break
    }
  }
  ppv = (rocRes$TP * prevalence[i])/(rocRes$TP * prevalence[i] + (1-prevalence[i]) * rocRes$FP)
  npv = ((1-prevalence[i]) * (1-rocRes$FP))/((1-prevalence[i]) * (1-rocRes$FP) + prevalence[i] * (1-rocRes$TP))
  rocRes$thrs[which.max(ppv + npv - 1)]
})

#The weibull scale and shape matters
minVisits = 5
patientDsList = split(testDs, testDs$P_ID)

res = vector("list", length(patientDsList))
tStart = Sys.time()
#for(i in 2:2){

res = foreach(i=201:250,
         .packages = c("splines", "JMbayes")) %dopar%{

  patientDs_i = patientDsList[[i]]
  
  patientDs_i$expectedFailureTime = NA
  patientDs_i$survTimeYouden = NA
  patientDs_i$survTimeMarkedNess = NA
  
  biopsy_times = list(expectedFailureTime = c(), survTimeYouden = c(), survTimeMarkedNess = c())
  
  for(j in minVisits:(timesPerSubject-1)){
    persTestDs = patientDs_i[patientDs_i$visitNumber <= j,]
    temp_lasttime = max(persTestDs$visitTimeYears)
    nearest_time_index = which(abs(dynamicCutOffTimes-temp_lasttime)==min(abs(dynamicCutOffTimes-temp_lasttime)))
    nearest_time_index = nearest_time_index[1]

    patientDs_i$expectedFailureTime[j] = expectedCondFailureTime(persTestDs)
    patientDs_i$survTimeYouden[j] = pDynSurvTime(survProb = youdenDynamicCutoffValues[nearest_time_index], persTestDs)
    patientDs_i$survTimeMarkedNess[j] = pDynSurvTime(survProb = markedNess[nearest_time_index], persTestDs)
    
    if(patientDs_i$expectedFailureTime[j] < patientDs_i$visitTimeYears[j + 1]){
      biopsy_times$expectedFailureTime = c(biopsy_times$expectedFailureTime, patientDs_i$expectedFailureTime[j])
    }
    
    if(patientDs_i$survTimeYouden[j] < patientDs_i$visitTimeYears[j + 1]){
      biopsy_times$survTimeYouden = c(biopsy_times$survTimeYouden, patientDs_i$survTimeYouden[j])
    }
    
    if(patientDs_i$survTimeMarkedNess[j] < patientDs_i$visitTimeYears[j + 1]){
      biopsy_times$survTimeMarkedNess = c(biopsy_times$survTimeMarkedNess, patientDs_i$survTimeMarkedNess[j])
    }
    print(j)
  }
  
  if(i > read.csv("Rdata/Gleason as event/simStudy.csv")[,2]){
    write.csv(i, file = "Rdata/Gleason as event/simStudy.csv") 
  }
  
  #res[[i]]=list(patientDs_i, biopsy_times)
  list(patientDs_i, biopsy_times)
}

tEnd = Sys.time()

save.image("Rdata/Gleason as event/simStudy1000.Rdata.RData")

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
