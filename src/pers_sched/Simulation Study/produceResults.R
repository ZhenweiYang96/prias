plotBiopsy2DPlot = function(nbSummary, offsetSummary, nbMeasure = "Median", offsetMeasure = "Median"){
  methodNames = colnames(nbSummary)
  methodCategory = sapply(methodNames, substr, 1, 10)
  
  df = data.frame(nb = nbSummary[nbMeasure,], offset=offsetSummary[offsetMeasure,], 
                  method=colnames(nbSummary), methodCategory = methodCategory)
  p = ggplot(data=df, aes(x=nb, y=offset*12, label=method)) + 
    geom_label(aes(fill = methodCategory), colour = "white", fontface = "bold") + 
    xlab(paste(nbMeasure,"Number of biopsies")) + ylab(paste(offsetMeasure, "Offset (months)")) + 
    scale_y_continuous(breaks=seq(0, 100, 1)) + scale_x_continuous(breaks=seq(0, 100, 0.5))
  print(p)
}

plot2DBoxPlot = function(nbSummary, offsetSummary){
  methodNames = colnames(nbSummary)
  methodCategory = sapply(methodNames, substr, 1, 10)
  
  df = data.frame(nb_x1 = nbSummary["1st Qu.",], nb_x2 = nbSummary["3rd Qu.",],
                  offset_y1 = offsetSummary["1st Qu.",], offset_y2 = offsetSummary["3rd Qu.",],
                  method=colnames(nbSummary), methodCategory = methodCategory)
  
  p = ggplot(data=df) + 
    geom_rect(mapping=aes(xmin=nb_x1, xmax=nb_x2, ymin=offset_y1*12, ymax=offset_y2*12, fill=methodCategory), color="black", alpha=0.5) +
    xlab("Number of biopsies") + ylab("Offset (months)")
  
  print(p)
}

produceResultImages = function(rDataFolder, simNumbers, DtSubFolder = "Dt_1", imgWidth=640, imgHeight=480){
  
  ct = makeCluster(detectCores())
  registerDoParallel(ct)
  
  patientCounts = foreach(simNum = simNumbers, .packages = c("ggplot2"),
                          .export=c("timesPerSubject", "plotBiopsy2DPlot"))%dopar%{
                            simFileLocation = paste("Rdata/Gleason as event/Sim Study/",rDataFolder, "/", DtSubFolder, "/simDs", simNum, ".Rdata", sep="")
                            load(simFileLocation)
                            print(paste("SimNum:", simNum, "loaded"))
                            
                            simulatedDsList = temp
                            
                            minVisitsMethod = paste("_minVisit_", 5, sep="")
                            
                            #DtSubFolder = paste(DtSubFolder, "sc_8")
                            #biopsyResults = data.frame(simulatedDsList[[1]]$biopsyTimes[simulatedDsList[[1]]$biopsyTimes$weibullScale==8,])
                            biopsyResults = data.frame(simulatedDsList[[1]]$biopsyTimes)
                            biopsyResults$nb = as.numeric(as.character(biopsyResults$nb))
                            biopsyResults$offset = as.numeric(as.character(biopsyResults$offset))
                            
                            biopsyResults$methodCategory = sapply(biopsyResults$methodName, substr, 1, 10)
                            
                            imageFolderPath = paste("images/sim study/", rDataFolder, "/", simulatedDsList[[1]]$seed,"/", sep = "")
                            boxplotOffsetPath = paste(imageFolderPath, "boxplot_offset/", sep="")
                            boxplotNbPath = paste(imageFolderPath, "boxplot_nb/", sep="")
                            nbVsOffsetMedianPath = paste(imageFolderPath, "nbVsOffset_Median/", sep="")
                            nbVsOffsetMeanPath = paste(imageFolderPath, "nbVsOffset_Mean/", sep="")
                            
                            dir.create(boxplotOffsetPath, recursive = T)
                            dir.create(boxplotNbPath, recursive = T)
                            dir.create(nbVsOffsetMedianPath, recursive = T)
                            dir.create(nbVsOffsetMeanPath, recursive = T)
                            
                            nbSummary = do.call(cbind, by(data = biopsyResults[, "nb"], INDICES = biopsyResults$methodName, summary))
                            offsetSummary = do.call(cbind, by(data = biopsyResults[, "offset"], INDICES = biopsyResults$methodName, summary))
                            
                            #Progression time
                            png(width=640, height=480, filename = paste(imageFolderPath, "progression_hist.png", sep=""))
                            p = ggplot(data=simulatedDsList[[1]]$testDs.id) + geom_histogram(aes(progression_time)) +
                              xlab("Progression Time (years)")
                            print(p)
                            dev.off()
                            
                            #boxplot offset
                            png(width=imgWidth, height=imgHeight, filename = paste(boxplotOffsetPath, DtSubFolder, minVisitsMethod,".png", sep=""))
                            p = ggplot(data = biopsyResults) +
                              geom_boxplot(aes( reorder(methodName, offset, FUN=median), offset*12, fill=methodCategory)) +
                              scale_y_continuous(breaks = 12*seq(0,20, by = 0.5)) +
                              ylab("Biopsy offset (months)") + xlab("Method")
                            print(p)
                            dev.off()
                            
                            #boxplot nb
                            png(width=imgWidth, height=imgHeight, filename = paste(boxplotNbPath, DtSubFolder, minVisitsMethod, ".png", sep=""))
                            p = ggplot(data = biopsyResults) +
                              geom_boxplot(aes( reorder(methodName, nb, FUN=median), nb, fill=methodCategory)) +
                              scale_y_continuous(breaks = seq(0,20, by = 1)) +
                              ylab("Number of biopsies") + xlab("Method")
                            print(p)
                            dev.off()
                            
                            #nbVsOffset Median
                            png(width=imgWidth, height=imgHeight, filename = paste(nbVsOffsetMedianPath, DtSubFolder, minVisitsMethod, ".png", sep=""))
                            plotBiopsy2DPlot(nbSummary, offsetSummary)
                            dev.off()
                            
                            #nbVsOffset Median
                            png(width=imgWidth, height=imgHeight, filename = paste(nbVsOffsetMeanPath, DtSubFolder, minVisitsMethod, ".png", sep=""))
                            plotBiopsy2DPlot(nbSummary, offsetSummary, "Mean", "Mean")
                            dev.off()
                            
                            return(c(nrow(biopsyResults)))
                          }
  
  print(patientCounts)
  
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
                      boxplotOffsetPath = paste(imageFolderPath, "boxplot_offset/", sep="")
                      boxplotNbPath = paste(imageFolderPath, "boxplot_nb/", sep="")
                      nbVsOffsetMedianPath = paste(imageFolderPath, "nbVsOffset_Median/", sep="")
                      nbVsOffsetMeanPath = paste(imageFolderPath, "nbVsOffset_Mean/", sep="")
                      
                      file.remove(paste(imageFolderPath, "progression_hist.png", sep=""))
                      file.remove(paste(boxplotOffsetPath, DtSubFolder, minVisitsMethod,".png", sep=""))
                      file.remove(paste(boxplotNbPath, DtSubFolder, minVisitsMethod,".png", sep=""))
                      file.remove(paste(nbVsOffsetMedianPath, DtSubFolder, minVisitsMethod,".png", sep=""))
                      file.remove(paste(nbVsOffsetMeanPath, DtSubFolder, minVisitsMethod,".png", sep=""))
                    }
  
  stopCluster(ct)
}