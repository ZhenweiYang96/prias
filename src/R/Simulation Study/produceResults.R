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