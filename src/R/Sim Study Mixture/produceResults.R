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

plotBiopsy2DPlot = function(nbSummary, offsetSummary, nbMeasure = "Median", offsetMeasure = "Median"){
  df = data.frame(nb = nbSummary[nbMeasure,], offset=offsetSummary[offsetMeasure,], 
                  methodCategory=colnames(nbSummary))
  p = ggplot(data=df, aes(x=nb, y=offset*12, label=methodCategory)) + 
    geom_label() + 
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
                            setwd("C:\\Users\\838035\\Google Drive\\PhD\\src\\prias")
                            simFileLocation = paste("Rdata/Gleason as event/Sim Study/",rDataFolder, "/", DtSubFolder, "/simDs", simNum, ".Rdata", sep="")
                            load(simFileLocation)
                            print(paste("SimNum:", simNum, "loaded"))
                            
                            simulatedDsList = temp
                            
                            minVisitsMethod = paste("_minVisit_", 5, sep="")
                            
                            #DtSubFolder = paste(DtSubFolder, "sc_4")
                            #biopsyResults = data.frame(simulatedDsList[[1]]$biopsyTimes[simulatedDsList[[1]]$biopsyTimes$weibullScale==4,])
                            biopsyResults = data.frame(simulatedDsList[[1]]$biopsyTimes)
                            biopsyResults$nb = as.numeric(as.character(biopsyResults$nb))
                            biopsyResults$offset = as.numeric(as.character(biopsyResults$offset))
                            
                            biopsyResults$methodCategory = sapply(biopsyResults$methodName, function(x){
                              if(x=="expectedFailureTime"){
                                return("Mean")
                              }else if(x=="medianFailureTime"){
                                return("Median")
                              }else if(x=="PRIAS"){
                                return("PRIAS")
                              }else if(x=="JH"){
                                return("Annual")
                              }else if(x=="youden"){
                                return("Youden")
                              }else if(x=="f1score"){
                                return("F1score")
                              }else{
                                return(as.character(x))
                              }
                            })
                            
                            imageFolderPath = paste("images/sim study/", rDataFolder, "/", simulatedDsList[[1]]$seed,"/", sep = "")
                            boxplotOffsetPath = paste(imageFolderPath, "boxplot_offset/", sep="")
                            boxplotNbPath = paste(imageFolderPath, "boxplot_nb/", sep="")
                            nbVsOffsetMedianPath = paste(imageFolderPath, "nbVsOffset_Median/", sep="")
                            nbVsOffsetMeanPath = paste(imageFolderPath, "nbVsOffset_Mean/", sep="")
                            
                            dir.create(boxplotOffsetPath, recursive = T)
                            dir.create(boxplotNbPath, recursive = T)
                            dir.create(nbVsOffsetMedianPath, recursive = T)
                            dir.create(nbVsOffsetMeanPath, recursive = T)
                            
                            nbSummary = do.call(cbind, by(data = biopsyResults[, "nb"], INDICES = biopsyResults$methodCategory, summary))
                            offsetSummary = do.call(cbind, by(data = biopsyResults[, "offset"], INDICES = biopsyResults$methodCategory, summary))
                            
                            #Progression time
                            png(width=640, height=480, filename = paste(imageFolderPath, "progression_hist.png", sep=""))
                            p = ggplot(data=simulatedDsList[[1]]$testDs.id) + geom_histogram(aes(progression_time)) +
                              xlab("Progression Time (years)")
                            print(p)
                            dev.off()
                            
                            #boxplot offset
                            png(width=imgWidth, height=imgHeight, filename = paste(boxplotOffsetPath, DtSubFolder, minVisitsMethod,".png", sep=""))
                            p = ggplot(data = biopsyResults) +
                              geom_boxplot(aes( reorder(methodCategory, offset, FUN=mean), offset*12)) +
                              scale_y_continuous(breaks = 12*seq(0,20, by = 1)) +
                              ylab("Biopsy offset (months)") + xlab("Method")
                            print(p)
                            dev.off()
                            
                            #boxplot nb
                            png(width=imgWidth, height=imgHeight, filename = paste(boxplotNbPath, DtSubFolder, minVisitsMethod, ".png", sep=""))
                            p = ggplot(data = biopsyResults) +
                              geom_boxplot(aes( reorder(methodCategory, nb, FUN=mean), nb)) +
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

boxplotAllPatients = function(rDataFolder, simNumbers, DtSubFolder = "Dt_1", subPopulationWeibullScale = NA, subpopName = "all"){
  if(is.na(subPopulationWeibullScale)){
    stop("enter sub population weibull scale, or vector of scales")
  }
  
  ct = makeCluster(detectCores())
  registerDoParallel(ct)
  
  biopsyResults = foreach(simNum = simNumbers, .packages = c("ggplot2"),
                           .export=c("timesPerSubject"), .combine="rbind")%dopar%{
                             
                             setwd("C:\\Users\\838035\\Google Drive\\PhD\\src\\prias")
                             simFileLocation = paste("Rdata/Gleason as event/Sim Study/",rDataFolder, "/", DtSubFolder, "/simDs", simNum, ".Rdata", sep="")
                             load(simFileLocation)
                             
                             simulatedDsList = temp
                             
                             simulatedDsList[[1]]$biopsyTimes = simulatedDsList[[1]]$biopsyTimes[simulatedDsList[[1]]$biopsyTimes$weibullScale %in% subPopulationWeibullScale,]
                             simulatedDsList[[1]]$biopsyTimes[,"methodCategory"] = sapply(simulatedDsList[[1]]$biopsyTimes[,"methodName"], function(x){
                               if(x=="expectedFailureTime"){
                                 return("Expec. Time GR")
                               } 
                               
                               if(x=="medianFailureTime"){
                                 return("Median Time GR")
                               } 
                               
                               if(x=="PRIAS"){
                                 return("PRIAS")
                               }
                               
                               if(x=="JH"){
                                 return("Annual")
                               }
                               
                               if(x=="youden"){
                                 return("Youden's J")
                               }
                               
                               if(x=="f1score"){
                                 return("F1-Score")
                               }
                               
                               if(x=="MixedYouden"){
                                 return("Mixed Approach")
                               }
                             })
                             return(simulatedDsList[[1]]$biopsyTimes)
                           }
  
  biopsyResults = data.frame(biopsyResults)
  biopsyResults$nb = as.numeric(as.character(biopsyResults$nb))
  biopsyResults$offset = as.numeric(as.character(biopsyResults$offset)) * 12
  
  
  png(width=640, height=480, filename = paste("report/pers_schedule/images/sim_study/", "nbBoxPlot_",subpopName,".png", sep=""))
  p = ggplot(data = biopsyResults) +
    geom_boxplot(aes(reorder(methodCategory, nb, FUN=mean), nb)) +
    scale_y_continuous(breaks = seq(0,20, by = 1)) +
    ylab("Number of biopsies") + xlab("Method")
  print(p)
  dev.off()
  
  png(width=640, height=480, filename = paste("report/pers_schedule/images/sim_study/", "offsetBoxPlot_",subpopName,".png", sep=""))
  p = ggplot(data = biopsyResults) +
    geom_boxplot(aes( reorder(methodCategory, nb, FUN=mean), offset)) +
    scale_y_continuous(breaks = seq(0,240, by = 12)) +
    ylab("Biopsy offset (months)") + xlab("Method")
  print(p)
  dev.off()
  
  stopCluster(ct)
  
  return(biopsyResults)
}

poolInformation = function(rDataFolder, simNumbers, DtSubFolder = "Dt_1", subPopulationWeibullScale = NA, subpopName = "all"){
  if(is.na(subPopulationWeibullScale)){
    stop("enter sub population weibull scale, or vector of scales")
  }
  
  ct = makeCluster(detectCores())
  registerDoParallel(ct)
  
  resultsSummary = foreach(simNum = simNumbers, .packages = c("ggplot2"),
                          .export=c("timesPerSubject"))%dopar%{
                            
                            setwd("C:\\Users\\838035\\Google Drive\\PhD\\src\\prias")
                            simFileLocation = paste("Rdata/Gleason as event/Sim Study/",rDataFolder, "/", DtSubFolder, "/simDs", simNum, ".Rdata", sep="")
                            load(simFileLocation)
                            
                            simulatedDsList = temp
                            
                            #biopsyResults = data.frame(simulatedDsList[[1]]$biopsyTimes)
                            biopsyResults = data.frame(simulatedDsList[[1]]$biopsyTimes[simulatedDsList[[1]]$biopsyTimes$weibullScale %in% subPopulationWeibullScale,])
                            biopsyResults$nb = as.numeric(as.character(biopsyResults$nb))
                            biopsyResults$offset = as.numeric(as.character(biopsyResults$offset))
                            
                            biopsyResults$methodCategory = sapply(biopsyResults$methodName, function(x){
                              if(x=="expectedFailureTime"){
                                return("Expec. Time GR")
                              } 
                              
                              if(x=="medianFailureTime"){
                                return("Median Time GR")
                              } 
                              
                              if(x=="PRIAS"){
                                return("PRIAS")
                              }
                              
                              if(x=="JH"){
                                return("Annual")
                              }
                              
                              if(x=="youden"){
                                return("Youden's J")
                              }
                              
                              if(x=="f1score"){
                                return("F1-Score")
                              }
                              
                              if(x=="MixedYouden"){
                                return("Mixed Approach")
                              }
                            })
                            
                            #Now we need total eligible patients per method, nb Mean, Var, and Offset mean, var per method
                            totalPatientsPerMethod = by(biopsyResults$nb, biopsyResults$methodCategory, length)
                            nbMeanPerMethod = by(biopsyResults$nb, biopsyResults$methodCategory, mean)
                            nbVarPerMethod = by(biopsyResults$nb, biopsyResults$methodCategory, var)
                            offsetMeanPerMethod = by(biopsyResults$offset * 12, biopsyResults$methodCategory, mean)
                            offsetVarPerMethod = by(biopsyResults$offset * 12, biopsyResults$methodCategory, var)
                            
                            return(list(totalPatientsPerMethod = totalPatientsPerMethod, 
                                        nbMeanPerMethod = nbMeanPerMethod,
                                        nbVarPerMethod = nbVarPerMethod,
                                        offsetMeanPerMethod = offsetMeanPerMethod,
                                        offsetVarPerMethod = offsetVarPerMethod))
                            
                          }
  stopCluster(ct)
  
  methodNames = names(resultsSummary[[1]][[1]])
  paramNames = names(resultsSummary[[1]])
  
  finalResultSummary = matrix(data = NA, nrow = length(methodNames), ncol=length(paramNames))
  rownames(finalResultSummary) = methodNames
  colnames(finalResultSummary) = paramNames
  
  finalResultSummary[,"totalPatientsPerMethod"] = apply(sapply(resultsSummary, FUN = function(x){x[["totalPatientsPerMethod"]]}), MARGIN = 1, sum)
  finalResultSummary[,"nbMeanPerMethod"] = apply(sapply(resultsSummary, FUN = function(x){x[["nbMeanPerMethod"]] * x[["totalPatientsPerMethod"]]}), MARGIN = 1, FUN = sum) / finalResultSummary[,"totalPatientsPerMethod"]
  finalResultSummary[,"offsetMeanPerMethod"] = apply(sapply(resultsSummary, FUN = function(x){x[["offsetMeanPerMethod"]] * x[["totalPatientsPerMethod"]]}), MARGIN = 1, FUN = sum) / finalResultSummary[,"totalPatientsPerMethod"]
  finalResultSummary[,"nbVarPerMethod"] = apply(sapply(resultsSummary, FUN = function(x){x[["nbVarPerMethod"]] * (x[["totalPatientsPerMethod"]]-1)}), MARGIN = 1, FUN = sum) / (finalResultSummary[,"totalPatientsPerMethod"] - length(resultsSummary))
  finalResultSummary[,"offsetVarPerMethod"] = apply(sapply(resultsSummary, FUN = function(x){x[["offsetVarPerMethod"]] * (x[["totalPatientsPerMethod"]]-1)}), MARGIN = 1, FUN = sum) / (finalResultSummary[,"totalPatientsPerMethod"] - length(resultsSummary))

  png(width=640, height=480, filename = paste("report/pers_schedule/images/sim_study/", "meanNbVsOffset_",subpopName,".png", sep=""))
  p = qplot(x = finalResultSummary[,"nbMeanPerMethod"], y=finalResultSummary[,"offsetMeanPerMethod"], label=rownames(finalResultSummary), geom="label",
            xlab="Mean number of biopsies", ylab="Mean offset (months)", xlim=c(min(finalResultSummary[,"nbMeanPerMethod"])-0.5,max(finalResultSummary[,"nbMeanPerMethod"])+0.25))
  print(p)
  dev.off()
  
  write.csv(finalResultSummary, file = paste("report/pers_schedule/csv/", "pooledInfo_", subpopName ,".csv", sep=""))
  
  return(finalResultSummary)
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