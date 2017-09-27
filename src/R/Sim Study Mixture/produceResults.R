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

boxplotAllPatients = function(rDataFolder, simNumbers, DtSubFolder = "Dt_1", subPopulationWeibullScale = NULL, subpopName = "all"){
  if(is.null(subPopulationWeibullScale)){
    stop("enter sub population weibull scale, or vector of scales")
  }
  
  ct = makeCluster(detectCores())
  registerDoParallel(ct)
  
  biopsyResults = foreach(simNum = simNumbers, .packages = c("ggplot2"),
                           .export=c("timesPerSubject"), .combine="rbind")%dopar%{
                             
                             setwd("C:\\Users\\838035\\Old Sim Results")
                             simFileLocation = paste("Rdata/Gleason as event/Sim Study/",rDataFolder, "/", DtSubFolder, "/simDs", simNum, ".Rdata", sep="")
                             load(simFileLocation)
                             
                             return(temp[[1]]$biopsyTimes)
                  }
  stopCluster(ct)
  
  png(width=640, height=480, filename = paste("report/pers_schedule/biometrics_submission/images/sim_study/", "nbBoxPlot_",subpopName,".png", sep=""))
  p = ggplot(data = biopsyResults[!biopsyResults$methodName %in% c("Youden", "Hybrid-Youden"),]) +
    geom_boxplot(aes(reorder(methodName, nb, FUN=mean), nb)) +
    scale_y_continuous(breaks = seq(0,20, by = 2)) +
    ylab("Number of biopsies") + xlab("Schedule") + 
    theme(text = element_text(size=14), axis.text=element_text(size=13))+ coord_flip()
  print(p)
  dev.off()
  
  png(width=640, height=480, filename = paste("report/pers_schedule/images/sim_study/", "offsetBoxPlot_",subpopName,".png", sep=""))
  p = ggplot(data = biopsyResults[!biopsyResults$methodName %in% c("Youden", "Hybrid-Youden"),]) +
    geom_boxplot(aes(reorder(methodName, nb, FUN=mean), offset)) +
    ylab("Biopsy offset (months)") + xlab("Schedule") + 
    theme(text = element_text(size=14), axis.text=element_text(size=13))+ coord_flip()
  print(p)
  dev.off()
  
  return(biopsyResults)
}

poolInformation = function(rDataFolder, simNumbers, DtSubFolder = "Dt_1", subPopulationWeibullScale = NULL, subpopName = "all"){
  if(is.null(subPopulationWeibullScale)){
    stop("enter sub population weibull scale, or vector of scales")
  }
  
  ct = makeCluster(detectCores())
  registerDoParallel(ct)
  
  resultsSummary = foreach(simNum = simNumbers, .packages = c("ggplot2"),
                          .export=c("timesPerSubject"))%dopar%{
                            
                            setwd("C:\\Users\\838035\\Old Sim Results")
                            simFileLocation = paste("Rdata/Gleason as event/Sim Study/",rDataFolder, "/", DtSubFolder, "/simDs", simNum, ".Rdata", sep="")
                            load(simFileLocation)
                            
                            biopsyResults = temp[[1]]$biopsyTimes[temp[[1]]$biopsyTimes$weibullScale %in% subPopulationWeibullScale,]
                            
                            #Now we need total eligible patients per method, nb Mean, Var, and Offset mean, var per method
                            totalPatientsPerMethod = by(biopsyResults$nb, biopsyResults$methodName, length)
                            nbMeanPerMethod = by(biopsyResults$nb, biopsyResults$methodName, mean)
                            nbVarPerMethod = by(biopsyResults$nb, biopsyResults$methodName, var)
                            offsetMeanPerMethod = by(biopsyResults$offset, biopsyResults$methodName, mean)
                            offsetVarPerMethod = by(biopsyResults$offset, biopsyResults$methodName, var)
                            offsetNbCov = by(biopsyResults[, c("offset", "nb")], biopsyResults$methodName, var)
                            
                            return(list(totalPatientsPerMethod = totalPatientsPerMethod, 
                                        nbMeanPerMethod = nbMeanPerMethod,
                                        nbVarPerMethod = nbVarPerMethod,
                                        offsetMeanPerMethod = offsetMeanPerMethod,
                                        offsetVarPerMethod = offsetVarPerMethod, 
                                        offsetNbCov = offsetNbCov))
                            
                          }
  stopCluster(ct)
  
  methodNames = names(resultsSummary[[1]][[1]])
  # paramNames = c("totalPatientsPerMethod","nbMeanPerMethod","nbMeanSDPerMethod",
  #                "offsetMeanPerMethod","offsetMeanSDPerMethod","nbSDPerMethod",
  #                "nbVarSDPerMethod","offsetSDPerMethod","offsetVarSDPerMethod")
  paramNames = c("totalPatientsPerMethod","nbMeanPerMethod", "offsetMeanPerMethod",
                  "nbSDPerMethod", "offsetSDPerMethod")
  
  finalResultSummary = matrix(data = NA, nrow = length(methodNames), ncol=length(paramNames))
  rownames(finalResultSummary) = methodNames
  colnames(finalResultSummary) = paramNames
  
  # finalResultSummary[,paramNames[1]] = apply(sapply(resultsSummary, FUN = function(x){x[["totalPatientsPerMethod"]]}), MARGIN = 1, sum)
  # 
  # finalResultSummary[,paramNames[2]] = apply(sapply(resultsSummary, FUN = function(x){x[["nbMeanPerMethod"]] * x[["totalPatientsPerMethod"]]}), MARGIN = 1, FUN = sum) / finalResultSummary[,"totalPatientsPerMethod"]
  # finalResultSummary[,paramNames[3]] = apply(sapply(resultsSummary, FUN = function(x){x[["nbMeanPerMethod"]]}), MARGIN = 1, FUN = var)
  # 
  # finalResultSummary[,paramNames[4]] = apply(sapply(resultsSummary, FUN = function(x){x[["offsetMeanPerMethod"]] * x[["totalPatientsPerMethod"]]}), MARGIN = 1, FUN = sum) / finalResultSummary[,"totalPatientsPerMethod"]
  # finalResultSummary[,paramNames[5]] = apply(sapply(resultsSummary, FUN = function(x){x[["offsetMeanPerMethod"]]}), MARGIN = 1, FUN = var)
  # 
  # finalResultSummary[,paramNames[6]] = apply(sapply(resultsSummary, FUN = function(x){x[["nbVarPerMethod"]] * (x[["totalPatientsPerMethod"]]-1)}), MARGIN = 1, FUN = sum) / (finalResultSummary[,"totalPatientsPerMethod"] - length(resultsSummary))
  # finalResultSummary[,paramNames[7]] = apply(sapply(resultsSummary, FUN = function(x){x[["nbVarPerMethod"]]}), MARGIN = 1, FUN = var)
  # 
  # finalResultSummary[,paramNames[8]] = apply(sapply(resultsSummary, FUN = function(x){x[["offsetVarPerMethod"]] * (x[["totalPatientsPerMethod"]]-1)}), MARGIN = 1, FUN = sum) / (finalResultSummary[,"totalPatientsPerMethod"] - length(resultsSummary))
  # finalResultSummary[,paramNames[9]] = apply(sapply(resultsSummary, FUN = function(x){x[["offsetVarPerMethod"]]}), MARGIN = 1, FUN = var)

  finalResultSummary[,paramNames[1]] = apply(sapply(resultsSummary, FUN = function(x){x[["totalPatientsPerMethod"]]}), MARGIN = 1, sum)
  
  finalResultSummary[,paramNames[2]] = apply(sapply(resultsSummary, FUN = function(x){x[["nbMeanPerMethod"]] * x[["totalPatientsPerMethod"]]}), MARGIN = 1, FUN = sum) / finalResultSummary[,"totalPatientsPerMethod"]
  finalResultSummary[,paramNames[3]] = apply(sapply(resultsSummary, FUN = function(x){x[["offsetMeanPerMethod"]] * x[["totalPatientsPerMethod"]]}), MARGIN = 1, FUN = sum) / finalResultSummary[,"totalPatientsPerMethod"]
  finalResultSummary[,paramNames[4]] = sqrt(apply(sapply(resultsSummary, FUN = function(x){x[["nbVarPerMethod"]] * (x[["totalPatientsPerMethod"]]-1)}), MARGIN = 1, FUN = sum) / (finalResultSummary[,"totalPatientsPerMethod"] - length(resultsSummary)))
  finalResultSummary[,paramNames[5]] = sqrt(apply(sapply(resultsSummary, FUN = function(x){x[["offsetVarPerMethod"]] * (x[["totalPatientsPerMethod"]]-1)}), MARGIN = 1, FUN = sum) / (finalResultSummary[,"totalPatientsPerMethod"] - length(resultsSummary)))

  finalResultSummary[,c(3,5)] = finalResultSummary[,c(3,5)] * 12 
  
  #Calculate mahalanobis distance
  tt = lapply(resultsSummary, FUN = function(resItem){Map('*', resItem[["offsetNbCov"]], resItem[["totalPatientsPerMethod"]]-1)})
  pooledCovariance = Map('/', lapply(methodNames, FUN=function(i){Reduce("+",lapply(tt, FUN=function(x){x[[i]]}))}), finalResultSummary[,"totalPatientsPerMethod"] - length(resultsSummary))
  rm(tt)
  mahal = rep(NA, 8)
  for(i in 1:8){mahal[i] = mahalanobis(c(0,1), center=finalResultSummary[i, c(3,2)], cov = pooledCovariance[[i]])}
  names(mahal) = methodNames
  #+ theme(text = element_text(size=15)) + coord_flip()
  #+ theme(text = element_text(size=15), axis.text.x = element_text(angle = 45, hjust = 1))
  
  png(width=640, height=480, filename = paste("report/pers_schedule/images/sim_study/", "nbMeanBoxPlot_",subpopName,".png", sep=""))
  ydata = c(sapply(resultsSummary, function(x){x[["nbMeanPerMethod"]]}))
  
  xdata = rep(methodNames,length(resultsSummary))
  indicesIgnore = xdata=="Hybrid-Youden" | xdata=="Youden"
  xdata = xdata[!indicesIgnore]
  ydata = ydata[!indicesIgnore]
  
  p = qplot(y=ydata,x=reorder(xdata, ydata, FUN=mean), geom = "boxplot", ylab="Mean number of biopsies", xlab="Schedule") + ticksY(0, 10, 0.5) +
    theme(text = element_text(size=12), axis.text=element_text(size=12)) + coord_flip()
  print(p)
  dev.off()
  
  png(width=640, height=480, filename = paste("report/pers_schedule/images/sim_study/", "offsetMeanBoxPlot_",subpopName,".png", sep=""))
  ydata = c(sapply(resultsSummary, function(x){x[["offsetMeanPerMethod"]]})) * 12
  xdata = rep(methodNames,length(resultsSummary))
  p = qplot(y=ydata,x=reorder(xdata, ydata, FUN=mean), geom = "boxplot", ylab="Mean biopsy offset (months)", xlab="Schedule") + ticksY(0, 100, 1) +
    theme(text = element_text(size=12), axis.text=element_text(size=12))+ coord_flip()
  print(p)
  dev.off()
  
  png(width=640, height=480, filename = paste("report/pers_schedule/images/sim_study/", "nbVarBoxPlot_",subpopName,".png", sep=""))
  ydata = c(sapply(resultsSummary, function(x){x[["nbVarPerMethod"]]}))
  xdata = rep(methodNames,length(resultsSummary))
  p = qplot(y=sqrt(ydata),x=reorder(xdata, ydata, FUN=mean), geom = "boxplot", ylab="Standard Deviation of number of biopsies", xlab="Schedule") + 
    ticksY(0, 10, 0.5) +
    theme(text = element_text(size=12), axis.text=element_text(size=12)) + coord_flip()
  print(p)
  dev.off()
  
  png(width=640, height=480, filename = paste("report/pers_schedule/images/sim_study/", "offsetVarBoxPlot_",subpopName,".png", sep=""))
  ydata = c(sapply(resultsSummary, function(x){x[["offsetVarPerMethod"]]})) * 12 * 12
  xdata = rep(methodNames,length(resultsSummary))
  p = qplot(y=sqrt(ydata),x=reorder(xdata, ydata, FUN=mean), geom = "boxplot", ylab="Standard Deviation of biopsy offset (months)", xlab="Schedule") + 
    ticksY(0, 50, 2.5) +
    theme(text = element_text(size=12), axis.text=element_text(size=12)) + coord_flip()
  print(p)
  dev.off()
  
  png(width=600, height=400, filename= paste("report/pers_schedule/biometrics_submission/images/sim_study/", "meanNbVsOffset_",subpopName,".png", sep=""))
  p = qplot(x = finalResultSummary[c(1:3, 5:6, 8), "nbMeanPerMethod"], 
            y=finalResultSummary[c(1:3, 5:6, 8), "offsetMeanPerMethod"], 
            label=rownames(finalResultSummary)[c(1:3, 5:6, 8)],
            xlab="Mean number of biopsies", ylab="Mean biopsy offset (months)", 
            xlim=c(min(finalResultSummary[c(1:3, 5:6, 8),"nbMeanPerMethod"])-0.5,
                   max(finalResultSummary[c(1:3, 5:6, 8),"nbMeanPerMethod"])+0.25)) + 
    theme(text = element_text(size=15), axis.text=element_text(size=15)) + geom_label(size=5)
            
  print(p)
  dev.off()
  
  write.csv(round(finalResultSummary,2), file = paste("report/pers_schedule/csv/", "pooledInfo_", subpopName ,".csv", sep=""))
  
  return(resultsSummary)
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

