prias_long$psadt = unlist(by(data = prias_long,INDICES = prias_long$P_ID, function(patientDs){
  
  totalRows = nrow(patientDs)
  result = rep(NA, totalRows)
  
  if(totalRows>=4){
    for(j in 4:totalRows){
      if(!is.na(patientDs$log2psa[j])){
        curVisitTime = patientDs$visitTimeYears[j]
        
        betaReg = lm(log2psa~visitTimeYears, data = patientDs[1:j, ])$coefficients[2]
        
        result[j] = 1/betaReg
      }
    }
  }
  
  return(result)
}))


#First you find time of biopsy per patient
temp = by(data=prias_long, INDICES = prias_long$P_ID, function(patientDs){
  
  print(patientDs$P_ID[1])
  
  fixedBiopsyTimes = c(0, 1, 4, 7, 10, 15, 20, 25)
  
  biopsyTimes = patientDs$visitTimeYears[!is.na(patientDs$gleason)]
  
  nearestBiopsyTime = rep(NA, length(biopsyTimes))
  nearestBiopsyTimeAbsOffset = rep(NA, length(biopsyTimes))
  
  for(i in 1:length(biopsyTimes)){
    tempindex = which(abs(fixedBiopsyTimes-biopsyTimes[i])==min(abs(fixedBiopsyTimes-biopsyTimes[i])))
    nearestBiopsyTime[i] = fixedBiopsyTimes[i]
    nearestBiopsyTimeAbsOffset[i] = abs(biopsyTimes[i] - nearestBiopsyTime[i])
  }
  
  maxAllowedFixDiff = 3.5/12
  
  possibleAnnualBiopsyTimes = biopsyTimes[nearestBiopsyTimeAbsOffset > maxAllowedFixDiff]
  
  count = 0
  for(possibleAnnualBiopsyTime in possibleAnnualBiopsyTimes){
    timeWindowLeft = possibleAnnualBiopsyTime - 1
    
    psadtInWindow = patientDs$psadt[patientDs$visitTimeYears <= possibleAnnualBiopsyTime & patientDs$visitTimeYears>=timeWindowLeft]
    
    if(any(psadtInWindow > 0 & psadtInWindow <=10, na.rm = T)){
      count = count + 1
    }
  }
  
  return(count)
})
