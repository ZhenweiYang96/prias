summarizeBiopsyResults =  function(dsId, patientRowNum, minVisits, biopsyIfLessThanTime=1, 
                                   methodName = "expectedFailureTime", biopsyEveryKYears = NA){
  trueProgressionTime = simulatedDsList[[dsId]]$testDs.id$progression_time[patientRowNum]
  
  patientId = simulatedDsList[[dsId]]$testDs.id$P_ID[patientRowNum]
  
  patientDs_i = simulatedDsList[[dsId]]$testDs[simulatedDsList[[dsId]]$testDs$P_ID == patientId,]
  
  visitTimeYears = patientDs_i$visitTimeYears
  biopsyTimes = patientDs_i[, methodName]
  
  nb = 0
  biopsyTimeOffset = NA
  
  lastBiopsyTime = 0
  
  for(j in minVisits:timesPerSubject){
    curVisitTime = visitTimeYears[j]
    predBiopsyTime = biopsyTimes[j]
    
    biopsy_gap_needed = (curVisitTime - lastBiopsyTime) < 1
    
    condition1 = !is.na(predBiopsyTime) & ((predBiopsyTime - curVisitTime) <= biopsyIfLessThanTime)
    condition2 = !is.na(biopsyEveryKYears) & !is.na(predBiopsyTime) & ((predBiopsyTime - lastBiopsyTime) >= biopsyEveryKYears)
    if(biopsy_gap_needed==FALSE & (condition1 | condition2)){
      biopsyTimesOfInterest = biopsyTimes[visitTimeYears >= curVisitTime & visitTimeYears<=predBiopsyTime]
      biopsyIndexOfInterest = which.min(biopsyTimesOfInterest) - 1 + j
      
      nb = nb + 1
      biopsyTimeOffset = biopsyTimes[biopsyIndexOfInterest] - trueProgressionTime
      lastBiopsyTime = biopsyTimes[biopsyIndexOfInterest]
      
      if(biopsyTimeOffset >= 0){
        break
      }
    }
  }
  
  return(c(patientRowNum = patientRowNum, nb = nb, biopsyTimeOffset = biopsyTimeOffset))
}

summarizeBiopsyResults =  function(dsId, patientRowNum, minVisits, biopsyIfLessThanTime=1, 
                                   methodName = "expectedFailureTime", biopsyEveryKYears = NA){
  trueProgressionTime = simulatedDsList[[dsId]]$testDs.id$progression_time[patientRowNum]
  
  patientId = simulatedDsList[[dsId]]$testDs.id$P_ID[patientRowNum]
  
  patientDs_i = simulatedDsList[[dsId]]$testDs[simulatedDsList[[dsId]]$testDs$P_ID == patientId,]
  
  visitTimeYears = patientDs_i$visitTimeYears
  biopsyTimes = patientDs_i[, methodName]
  
  nb = 0
  biopsyTimeOffset = NA
  
  lastBiopsyTime = 0
  proposedBiopsyTime = Inf
  enforcedBiopsyTime = Inf
  
  for(j in minVisits:timesPerSubject){
    curVisitTime = visitTimeYears[j]
    proposedBiopsyTime = biopsyTimes[j]
    
    #Account for the biopsy that needed to be done
    if(curVisitTime > enforcedBiopsyTime){
      lastBiopsyTime = enforcedBiopsyTime
      nb = nb + 1
      biopsyTimeOffset = enforcedBiopsyTime - trueProgressionTime
      enforcedBiopsyTime = Inf
    }
    
    if(!is.na(biopsyTimeOffset) & biopsyTimeOffset > 0){
      break
    }
    
    #If the K years have passed since last biopsy
    if(!is.na(biopsyEveryKYears)){
      if((curVisitTime - lastBiopsyTime) >= biopsyEveryKYears){
        enforcedBiopsyTime = curVisitTime
      }
    }
    
    #Now deal with the proposed biopsy time
    if(!is.na(proposedBiopsyTime)){
      if(proposedBiopsyTime < enforcedBiopsyTime){
        # if((proposedBiopsyTime - lastBiopsyTime) <= 1){
        #   enforcedBiopsyTime = lastBiopsyTime + 1
        # }else if(proposedBiopsyTime - curVisitTime <= 1){
        #   enforcedBiopsyTime = proposedBiopsyTime
        # }
        
        if((proposedBiopsyTime - lastBiopsyTime) > 1){
          enforcedBiopsyTime = proposedBiopsyTime
        }
      }
    }
  }
  
  if(enforcedBiopsyTime < Inf){
    nb = nb + 1
    biopsyTimeOffset = enforcedBiopsyTime - trueProgressionTime
  }
  
  return(c(patientRowNum = patientRowNum, nb = nb, biopsyTimeOffset = biopsyTimeOffset))
}

getBiopsyResults = function(dsId, biopsyIfLessThanTime = 1, biopsyEveryKYears = NA , minVisits, methodNames = c("expectedFailureTime", "survTime85",
                                                                                                                "survTimeYouden", "survTimeAccuracy", "survTimeMaxTPR", "survTimeF1Score")){
  subjectCount = nrow(simulatedDsList[[dsId]]$testDs.id)
  
  biopsyResults = data.frame(patientRowNum = numeric(), nb = numeric(), biopsyTimeOffset = numeric(), methodName=factor(levels = methodNames))
  for(methodName in methodNames){
    summary = sapply(1:subjectCount, summarizeBiopsyResults, dsId = dsId, methodName = methodName, minVisits = minVisits, biopsyIfLessThanTime=biopsyIfLessThanTime)  
    biopsyResults = rbind(biopsyResults, data.frame(t(summary), methodName = rep(methodName,ncol(summary))))
  }
  
  if(!any(is.na(biopsyEveryKYears))){
    for(K in biopsyEveryKYears){
      for(methodName in methodNames){
        extMethodName = paste(methodName, K, sep="")
        summary = sapply(1:subjectCount, summarizeBiopsyResults, dsId = dsId, methodName = methodName, minVisits = minVisits, biopsyIfLessThanTime=biopsyIfLessThanTime, biopsyEveryKYears = K)
        biopsyResults = rbind(biopsyResults, data.frame(t(summary), methodName = rep(extMethodName,ncol(summary))))
      }
    }
  }
  
  #Fixed schedule
  priasSummary = sapply(1:subjectCount, function(patientRowNum){
    progressionTime = simulatedDsList[[dsId]]$testDs.id$progression_time[patientRowNum]
    
    if(progressionTime<=1){
      return(c(patientRowNum = patientRowNum, nb = 1, biopsyTimeOffset = 1-progressionTime))
    }else if(progressionTime<=4){
      return(c(patientRowNum = patientRowNum, nb = 2, biopsyTimeOffset = 4-progressionTime))
    }else if(progressionTime<=7){
      return(c(patientRowNum = patientRowNum, nb = 3, biopsyTimeOffset = 7-progressionTime))
    }else if(progressionTime<=10){
      return(c(patientRowNum = patientRowNum, nb = 4, biopsyTimeOffset = 10-progressionTime))
    }else if(progressionTime<=15){
      return(c(patientRowNum = patientRowNum, nb = 5, biopsyTimeOffset = 15-progressionTime))
    }else if(progressionTime<=20){
      return(c(patientRowNum = patientRowNum, nb = 6, biopsyTimeOffset = 20-progressionTime))
    }
  }, simplify = T)
  biopsyResults = rbind(biopsyResults, data.frame(t(priasSummary), methodName = rep("PRIAS",ncol(priasSummary))))
  
  #Johns hopkins schedule
  johnsHopkinsSummary =  sapply(1:subjectCount, function(patientRowNum){
    progressionTime = simulatedDsList[[dsId]]$testDs.id$progression_time[patientRowNum]
    detectionTime = ceiling(progressionTime)
    return(c(patientRowNum = patientRowNum, nb = detectionTime-1, biopsyTimeOffset = detectionTime-progressionTime))
  })
  biopsyResults = rbind(biopsyResults, data.frame(t(johnsHopkinsSummary), methodName = rep("johnsSummary",ncol(johnsHopkinsSummary))))
  
  biopsyResults
}


produceResultImages = function(rDataFolder, simNumbers, DtSubFolder = "Dt_1", minVisits = 5, imgWidth=1280, imgHeight=960){
  
  ct = makeCluster(detectCores())
  registerDoParallel(ct)
  
  foreach(simNum = simNumbers, .packages = c("ggplot2"),
          .export=c("getBiopsyResults", "summarizeBiopsyResults", 
                    "timesPerSubject", "plotBiopsy2DPlot"))%dopar%{
    simFileLocation = paste("Rdata/Gleason as event/Sim Study/",rDataFolder, "/", DtSubFolder, "/simDs", simNum, ".Rdata", sep="")
    load(simFileLocation)
    print(paste("SimNum:", simNum, "loaded"))
    
    simulatedDsList = temp
  
    minVisitsMethod = paste("_minVisit_", minVisits, sep="")

    biopsyResults = getBiopsyResults(1, biopsyIfLessThanTime = 1, biopsyEveryKYears=c(2,3), 
                                     minVisits = minVisits)

    incompleteRowNum = unique(biopsyResults[is.na(biopsyResults$biopsyTimeOffset) | biopsyResults$biopsyTimeOffset < 0, ]$patientRowNum)
    biopsyResultsCC = biopsyResults[!(biopsyResults$patientRowNum %in% incompleteRowNum),]
    biopsyResultsCC$methodCategory = sapply(biopsyResultsCC$methodName, substr, 1, 10)
    biopsyResultsCC$methodTrimmed = sapply(biopsyResultsCC$methodName, function(x){
      if(as.character(x) == "PRIAS"){
        return(x)
      }
      
      if(as.character(x) == "johnsSummary"){
        return(factor("Johns"))
      }
      
      if(substr(as.character(x),1,3)=="exp"){
        x = as.character(x)
        return(factor(paste("E(T)-", substr(x,nchar(x)-1,nchar(x)), sep="")))
      }
      
      x = as.character(x)
      lastPart = substr(x, nchar(x)-6, nchar(x))
      factor(lastPart)
    })
    
    imageFolderPath = paste("images/sim study/", rDataFolder, "/", simulatedDsList[[1]]$seed,"/", sep = "")
    boxplotOffsetPath = paste(imageFolderPath, "boxplot_offset/")
    boxplotNbPath = paste(imageFolderPath, "boxplot_nb/")
    nbVsOffsetMedianPath = paste(imageFolderPath, "nbVsOffset_Median/")
    nbVsOffsetMeanPath = paste(imageFolderPath, "nbVsOffset_Mean/")
    
    dir.create(boxplotOffsetPath, recursive = T)
    dir.create(boxplotNbPath, recursive = T)
    dir.create(nbVsOffsetMedianPath, recursive = T)
    dir.create(nbVsOffsetMeanPath, recursive = T)
    
    nbSummary = do.call(cbind, by(data = biopsyResultsCC[, "nb"], INDICES = biopsyResultsCC$methodName, summary))
    offsetSummary = do.call(cbind, by(data = biopsyResultsCC[, "biopsyTimeOffset"], INDICES = biopsyResultsCC$methodName, summary))
    
    #Progression time
    png(width=640, height=480, filename = paste(imageFolderPath, "progression_hist.png", sep=""))
    p = ggplot(data=simulatedDsList[[1]]$testDs.id) + geom_histogram(aes(progression_time)) +
      xlab("Progression Time (years)")
    print(p)
    dev.off()
    
    #boxplot offset
    png(width=imgWidth, height=imgHeight, filename = paste(boxplotOffsetPath, DtSubFolder, minVisitsMethod,".png", sep=""))
    p = ggplot(data = biopsyResultsCC) +
      geom_boxplot(aes( reorder(methodTrimmed, biopsyTimeOffset, FUN=median), biopsyTimeOffset*12, fill=methodCategory)) +
      scale_y_continuous(breaks = 12*seq(0,20, by = 0.25)) +
      ylab("Biopsy offset (months)") + xlab("Method")
    print(p)
    dev.off()
    
    #boxplot nb
    png(width=imgWidth, height=imgHeight, filename = paste(boxplotNbPath, DtSubFolder, minVisitsMethod,".png", sep=""))
    p = ggplot(data = biopsyResultsCC) +
      geom_boxplot(aes( reorder(methodTrimmed, nb, FUN=median), nb, fill=methodCategory)) +
      scale_y_continuous(breaks = seq(0,20, by = 1)) +
      ylab("Number of biopsies") + xlab("Method")
    print(p)
    dev.off()
   
    #nbVsOffset Median
    png(width=imgWidth, height=imgHeight, filename = paste(nbVsOffsetMedianPath, DtSubFolder, minVisitsMethod,".png", sep=""))
    plotBiopsy2DPlot(nbSummary, offsetSummary)
    dev.off()
    
    #nbVsOffset Median
    png(width=imgWidth, height=imgHeight, filename = paste(nbVsOffsetMeanPath, DtSubFolder, minVisitsMethod,".png", sep=""))
    plotBiopsy2DPlot(nbSummary, offsetSummary, "Mean", "Mean")
    dev.off()
  }
  
  stopCluster(ct)
}


deleteResultImages = function(rDataFolder, simNumbers, DtSubFolder = "Dt_1", minVisits = 1, imgWidth=1280, imgHeight=960){
  
  ct = makeCluster(detectCores())
  registerDoParallel(ct)
  
  foreach(simNum = simNumbers, .packages = c("ggplot2"),
          .export=c("getBiopsyResults", "summarizeBiopsyResults", 
                    "timesPerSubject", "plotBiopsy2DPlot"))%dopar%{
                      simFileLocation = paste("Rdata/Gleason as event/Sim Study/",rDataFolder, "/", DtSubFolder, "/simDs", simNum, ".Rdata", sep="")
                      load(simFileLocation)
                      print(paste("SimNum:", simNum, "loaded"))
                      
                      simulatedDsList = temp
                      
                      minVisitsMethod = paste("_minVisit_", minVisits, sep="")
                      
                      imageFolderPath = paste("images/sim study/", rDataFolder, "/", simulatedDsList[[1]]$seed,"/", sep = "")
                      boxplotOffsetPath = paste(imageFolderPath, "boxplot_offset/")
                      boxplotNbPath = paste(imageFolderPath, "boxplot_nb/")
                      nbVsOffsetMedianPath = paste(imageFolderPath, "nbVsOffset_Median/")
                      nbVsOffsetMeanPath = paste(imageFolderPath, "nbVsOffset_Mean/")
                      
                      file.remove(paste(imageFolderPath, "progression_hist.png", sep=""))
                      file.remove(paste(boxplotOffsetPath, DtSubFolder, minVisitsMethod,".png", sep=""))
                      file.remove(paste(boxplotNbPath, DtSubFolder, minVisitsMethod,".png", sep=""))
                      file.remove(paste(nbVsOffsetMedianPath, DtSubFolder, minVisitsMethod,".png", sep=""))
                      file.remove(paste(nbVsOffsetMeanPath, DtSubFolder, minVisitsMethod,".png", sep=""))
                    }
  
  stopCluster(ct)
}