
switchToAnnual = function(ds){
  psaDt = 1/(lm(log2psa~visitTimeYears, data = ds)$coefficients[2])
  return(psaDt>=0 & psaDt<=10)
}

#psa_ds = prias_long[!is.na(prias_long$psa) & prias_long$P_ID %in% 132,]
psa_ds = prias_long[!is.na(prias_long$psa),]
#gleason_ds = prias_long[!is.na(prias_long$gleason) & prias_long$P_ID %in% 132,]
gleason_ds = prias_long[!is.na(prias_long$gleason),]

psa_ds$P_ID = droplevels(psa_ds$P_ID)
gleason_ds$P_ID = droplevels(gleason_ds$P_ID)

temp = by(psa_ds$P_ID, data=psa_ds, FUN = function(patientDs){
  print(patientDs$P_ID[1])
  lastBiopsyTime = 0
  proposedBiopsyTime = Inf
  
  proposedAndDoneList = list()
  
  fixedSchedule = c(1, 4, 7, 10, 15, 20)
  
  isAnnualBiopsyProposed = F
  proposedBiopsyDone = F
  
  gleason_ds_pat = gleason_ds[gleason_ds$P_ID == patientDs$P_ID[1],]
  usedIndices = c()
  
  minIndex = 1
  #Min number of measurements before which PSA-DT can't be used
  for(j in 1:nrow(patientDs)){
    curVisitTime = patientDs$visitTimeYears[j]
    
    if(proposedBiopsyDone==F | curVisitTime >= gleason_ds_pat$visitTimeYears[minIndex]){
      if(curVisitTime >= proposedBiopsyTime & proposedBiopsyDone==F & 
         any(sapply(proposedAndDoneList,FUN = function(x){
           r1 = x["proposedBiopsyTime"]==proposedBiopsyTime
           r2 = isAnnualBiopsyProposed & round(proposedBiopsyTime - proposedAndDoneList[[length(proposedAndDoneList)]]["proposedBiopsyTime"],2)<1
           return(r1 | r2)
         }))==F){
        
        if(isAnnualBiopsyProposed==F){
          minDiffUpperThres = 0.5
          minDiffLowerThres = -0.5
        }else{
          minDiffUpperThres = 1
          minDiffLowerThres = -0.5
        }
        
        minIndex = which.min(abs(gleason_ds_pat$visitTimeYears - proposedBiopsyTime))
        minDiff = (gleason_ds_pat$visitTimeYears - proposedBiopsyTime)[minIndex]
        if(minDiff > minDiffLowerThres & minDiff < minDiffUpperThres & !(minIndex %in% usedIndices)){
          usedIndices = c(usedIndices, minIndex)
          
          proposedAndDoneList[[length(proposedAndDoneList)+1]] = c("proposedBiopsyTime"=proposedBiopsyTime, 
                                                                   "conductedBiopsyTime" = gleason_ds_pat$visitTimeYears[minIndex],
                                                                   "isAnnualBiopsyProposed"=isAnnualBiopsyProposed)
          
          proposedBiopsyDone = T
          proposedBiopsyTime = Inf
          next
        }else{
          proposedAndDoneList[[length(proposedAndDoneList)+1]] = c("proposedBiopsyTime"=proposedBiopsyTime, 
                                                                   "conductedBiopsyTime" = NA,
                                                                   "isAnnualBiopsyProposed"=isAnnualBiopsyProposed)
          
          proposedBiopsyTime = Inf
        }
      }
      
      proposedBiopsyDone = F
      
      lastBiopsyTime = max(0,max(gleason_ds_pat$visitTimeYears[gleason_ds_pat$visitTimeYears<=curVisitTime]))
      
      if(j>=4 & switchToAnnual(patientDs[1:j, ])){
        # if((curVisitTime - lastBiopsyTime) >= 1){
        #   proposedBiopsyTime = curVisitTime
        # }else{
        #   proposedBiopsyTime = lastBiopsyTime + 1
        # }
        proposedBiopsyTime = curVisitTime
        isAnnualBiopsyProposed = T
      }else{
        proposedBiopsyTime = fixedSchedule[which(fixedSchedule >= curVisitTime)[1]]
        
        #if(proposedBiopsyTime - lastBiopsyTime < 1){
        #  proposedBiopsyTime = lastBiopsyTime + 1
        #}
        isAnnualBiopsyProposed = F
      }
    }
  }
  
  #check for 287 why it enters when he is not eligible for biopsy at year 10
  if(proposedBiopsyTime!=Inf & any(sapply(proposedAndDoneList,FUN = function(x){
    r1 = x["proposedBiopsyTime"]==proposedBiopsyTime
    r2 = isAnnualBiopsyProposed & round(proposedBiopsyTime - proposedAndDoneList[[length(proposedAndDoneList)]]["proposedBiopsyTime"],2)<1
    return(r1 | r2)
  }))==F){
    if(isAnnualBiopsyProposed==F){
      minDiffUpperThres = 0.5
      minDiffLowerThres = -0.5
    }else{
      minDiffUpperThres = 1
      minDiffLowerThres = -0.5
    }
    
    minIndex = which.min(abs(gleason_ds_pat$visitTimeYears - proposedBiopsyTime))
    minDiff = (gleason_ds_pat$visitTimeYears - proposedBiopsyTime)[minIndex]
    if(minDiff > minDiffLowerThres & minDiff < minDiffUpperThres & !(minIndex %in% usedIndices)){
      usedIndices = c(usedIndices, minIndex)
      
      proposedAndDoneList[[length(proposedAndDoneList)+1]] = c("proposedBiopsyTime"=proposedBiopsyTime, 
                                                               "conductedBiopsyTime" = gleason_ds_pat$visitTimeYears[minIndex],
                                                               "isAnnualBiopsyProposed"=isAnnualBiopsyProposed)
      
    }
  }
  
  print(length(usedIndices)==(nrow(gleason_ds_pat)-1))
  
  return(proposedAndDoneList)
})

year1=unlist(sapply(temp, function(x){
  retVec = c()
  if(length(x)>0){
    for(i in 1:length(x)){
      #if(x[[i]]["isAnnualBiopsyProposed"] & x[[i]]["proposedBiopsyTime"]>5.5){
      if(x[[i]]["isAnnualBiopsyProposed"] & x[[i]]["proposedBiopsyTime"]>8.5 & x[[i]]["proposedBiopsyTime"]<=9.5){
        retVec = c(retVec, x[[i]]["conductedBiopsyTime"])
      }
    }
    if(length(retVec)>0)
      return(retVec)
    else
      return(-1)
  }else{
    return(-1)
  }
}))
year1 = year1[!(year1 %in% c(-1))]

table(!is.na(year1))
table(!is.na(year1))/length(year1)

# randomPid = sample(prias.id$P_ID,size = 1)
# temp[[as.character(randomPid)]]
# gleason_ds$visitTimeYears[gleason_ds$P_ID==randomPid]
# psa_ds$visitTimeYears[psa_ds$P_ID==randomPid]
# 

###########################
## Now calculate in a simulation study how non compliance affects PRIAS schedule
###########################

library(doParallel)

resFiles = list.files("/home/a_tomer/Results/final_res_2nd_paper/", full.names = T)

ct = makeCluster(4)
registerDoParallel(ct)
prias_real_bokhorst = foreach(i=1:length(resFiles), .combine="rbind") %dopar%{
  
  load(resFiles[i])
  
  seed = jointModelData$seed
  set.seed(seed)
  
  runPRIAS_RealSchedule = function(testDs.id, testDs){
    fixedSchedule = c(1, 4, 7, 10)
    
    psa = 2^testDs$log2psaplus1 - 1
    psa[psa<=0] = 0.01
    testDs$log2psa = log(psa, base = 2)
    
    getComplianceRate = function(isAnnualBiopsy, biopsyTime){
      if(abs(biopsyTime - 1)<=0.5){
        return(0.8132875)
      }else if(abs(biopsyTime - 4)<=0.5){
        return(0.5787037)
      }else if(abs(biopsyTime - 7)<=0.5){
        return(0.5)
      }else if(abs(biopsyTime - 10)<=0.5){
        return(0.2727273)
      }else if(isAnnualBiopsy & biopsyTime > 1.5 & biopsyTime<=2.5){
        return(0.3755102)
      }else if(isAnnualBiopsy & biopsyTime > 2.5 & biopsyTime<=3.5){
        return(0.4911591)
      }else if(isAnnualBiopsy & biopsyTime > 4.5 & biopsyTime<=5.5){
        return(0.2989691)
      }else if(isAnnualBiopsy & biopsyTime > 5.5 & biopsyTime<=6.5){
        return(0.364486)
      }else if(isAnnualBiopsy & biopsyTime > 7.5 & biopsyTime<=8.5){
        return(0.06666667)
      }else if(isAnnualBiopsy & biopsyTime > 8.5 & biopsyTime<=9.5){
        return(0.25)
      }else{
        return(1)
      }
    }
    
    getComplianceRate_PRIAS = function(isAnnualBiopsy, biopsyTime){
      if(biopsyTime>0 & biopsyTime<=1){
        return(0.81)
      }else if(biopsyTime>1 & biopsyTime<=2){
        return(0.24)
      }else if(biopsyTime>2 & biopsyTime<=3){
        return(0.29)
      }else if(biopsyTime>3 & biopsyTime<=4){
        return(0.60)
      }else if(biopsyTime>4 & biopsyTime<=5){
        return(0.18)
      }else if(biopsyTime>5 & biopsyTime<=6){
        return(0.21)
      }else if(biopsyTime>6 & biopsyTime<=7){
        return(0.53)
      }else if(biopsyTime>7 & biopsyTime<=8){
        return(0.09)
      }else if(biopsyTime>8 & biopsyTime<=9){
        return(0.0)
      }else{
        return(0.33)
      }
    }
    
    switchToAnnual = function(ds){
      psaDt = 1/(lm(log2psa~visitTimeYears, data = ds)$coefficients[2])
      return(psaDt>=0 & psaDt<=10)
    }
    
    getNbOffset = function(pid){
      patientDs = testDs[testDs$P_ID == pid,]
      progressionTime = testDs.id$progression_time[testDs.id$P_ID == pid]
      
      nb = 0
      offset = NA
      
      #4 is the minimum number of measurements before which PSA-DT can't be used
      #The biopsy at year 1 will be able to detect the prostate cancer
      if(nrow(patientDs)<4){
        nb = 1
        offset = 1 - progressionTime
      }else{
        lastBiopsyTime = 0
        proposedBiopsyTime = Inf
        
        isAnnualSchedule = F
        
        #Min number of measurements before which PSA-DT can't be used
        for(j in 4:nrow(patientDs)){
          curVisitTime = patientDs$visitTimeYears[j]
          
          if(curVisitTime >= proposedBiopsyTime){
            complianceRate = getComplianceRate_PRIAS(isAnnualSchedule, proposedBiopsyTime)
            comply = rbinom(1,1,complianceRate)
            #print(comply)
            if(comply==T){
              nb = nb + 1
              offset = proposedBiopsyTime - progressionTime
              lastBiopsyTime = proposedBiopsyTime
              proposedBiopsyTime = Inf
            }
          }
          
          if(!is.na(offset) & offset > 0){
            break
          }
          
          if(switchToAnnual(patientDs[1:j, ])){
            if((curVisitTime - lastBiopsyTime) >= 1){
              proposedBiopsyTime = curVisitTime
            }else{
              proposedBiopsyTime = lastBiopsyTime + 1
            }
            
            isAnnualSchedule = T
          }else{
            proposedBiopsyTime = fixedSchedule[which(fixedSchedule >= curVisitTime)[1]]
            
            if(proposedBiopsyTime - lastBiopsyTime < 1){
              proposedBiopsyTime = lastBiopsyTime + 1
            }
            
            isAnnualSchedule = F
          }
        }
        
        #If we iterated through the entire vector of visit times
        if(proposedBiopsyTime < Inf){
          complianceRate = getComplianceRate(isAnnualSchedule, proposedBiopsyTime)
          comply = rbinom(1,1,complianceRate)
          #print(comply)
          if(comply==T){
            nb = nb + 1
            offset = proposedBiopsyTime - progressionTime
          }
        }
      }
      
      return(c(progressionTime, nb, offset))
    }
    
    return(t(sapply(testDs.id$P_ID, getNbOffset)))
  }
  
  prias_real_schedule =runPRIAS_RealSchedule(jointModelData$testData$testDs.id, jointModelData$testData$testDs)
  
  rm(jointModelData)
  #save(prias_real_schedule, file=paste("/home/a_tomer/Results/prias_real/prias_real_schedule_", seed, ".Rdata", sep = ""))
  return(prias_real_schedule)
}

stopCluster(ct)
