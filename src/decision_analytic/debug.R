auc1 = aucJM(mvJoint_dre_psa_dre_value, dre_psa_data_set, Tstart = 1, Thoriz=1 + 1.0101, idVar="P_ID")
save(auc1, file = "auc1.Rdata")
auc2 = aucJM(mvJoint_dre_psa_dre_value, dre_psa_data_set, Tstart = 2, Thoriz=2 + 1.0101, idVar="P_ID")
save(auc2, file = "auc2.Rdata")
auc3 = aucJM(mvJoint_dre_psa_dre_value, dre_psa_data_set, Tstart = 3, Thoriz=3 + 1.0101, idVar="P_ID")
save(auc3, file = "auc3.Rdata")


auc1 = aucJM(mvJoint_psa, psa_data_set, Tstart = 1, Thoriz=1 + 1.0101, idVar="P_ID")
save(auc1, file = "auc1_psa.Rdata")
auc2 = aucJM(mvJoint_psa, psa_data_set, Tstart = 2, Thoriz=2 + 1.0101, idVar="P_ID")
save(auc2, file = "auc2_psa.Rdata")
auc3 = aucJM(mvJoint_psa, psa_data_set, Tstart = 3, Thoriz=3 + 1.0101, idVar="P_ID")
save(auc3, file = "auc3_psa.Rdata")




library(doParallel)

resFiles = list.files("/home/a_tomer/Results/res20/", full.names = T)

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
