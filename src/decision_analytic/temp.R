
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
