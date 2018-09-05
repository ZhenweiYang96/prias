#load("Rdata/decision_analytic/DRE_PSA/mvJoint_dre_psa_dre_value.Rdata")
#load("Rdata/decision_analytic/cleandata.Rdata")
source("src/decision_analytic/load_lib.R")

getLastBiopsyTime = function(pid, lastnumber=1, upperLimitTime = Inf){
  temp = prias_long[prias_long$P_ID %in% pid & prias_long$visitTimeYears<=upperLimitTime,][, c("visitTimeYears", "gleason")]
  lastBiopsyTime = tail(temp[complete.cases(temp),]$visitTimeYears, lastnumber)[1]
  return(lastBiopsyTime)
}

plotDynamicRiskProbProfile = function(pid, fittedJointModel, maxVisitTime, 
                               maxPredictionTime=NA, lastBiopsyTime=NA, 
                               usecolor=F, threshold=0.15,
                               FONT_SIZE=12, POINT_SIZE = 2, DRE_PSA_Y_GAP=0.1,
                               LABEL_SIZE = 1){
  dataset = fittedJointModel$model_info$mvglmer_components$data
  
  patientDs = dataset[dataset$P_ID == pid & dataset$visitTimeYears<=maxVisitTime,]
  
  lastBiomarkerTime = max(patientDs$visitTimeYears)
  if(is.na(lastBiopsyTime)){
    lastBiopsyTime = getLastBiopsyTime(pid, upperLimitTime = maxVisitTime)
  }
  if(is.na(maxPredictionTime)){
    maxPredictionTime = lastBiomarkerTime
  }
  
  if(maxPredictionTime == lastBiopsyTime){
    stop("Max prediction time and last biopsy time are same")
  }
  
  futureTimes = seq(lastBiopsyTime, maxPredictionTime, length.out = 20)
  
  sfit = survfitJM(fittedJointModel, patientDs, idVar="P_ID", 
                   survTimes = futureTimes, last.time = lastBiopsyTime)
  
  patientDs$fitted_high_dre_prob = plogis(sfit$fitted.y[[1]]$high_dre)
  patientDs$fitted_log2psaplus1 = sfit$fitted.y[[1]]$log2psaplus1
  
  #The base of axes in this plot is of DRE
  minYLeft = 0
  maxYleft = 2 + DRE_PSA_Y_GAP * 2
  
  psaDs = patientDs[, c("visitTimeYears", "log2psaplus1", "fitted_log2psaplus1")]
  maxPSA = max(psaDs[,-1], na.rm=T)
  psaDs[,-1] = psaDs[,-1] / maxPSA + maxYleft/2 + DRE_PSA_Y_GAP
  
  riskDf = data.frame(sfit$summaries[[1]])
  riskDf = rbind(c("times"=futureTimes[1], "Mean"=1, "Median"=1, "Lower"=1, "Upper"=1), riskDf)
  # -1 because 1st column is time
  riskDf[,-1] = 1 - riskDf[,-1]
  maxMeanRiskLabel = paste0(round(max(riskDf$Mean)*100,2), "%")
  riskDf[,-1] = riskDf[, -1] * (maxYleft - minYLeft) + minYLeft
  maxMeanRisk = max(riskDf$Mean)
  
  xTicks = seq(0, maxPredictionTime, length.out = 6)
  xTicks = xTicks[abs(xTicks - lastBiopsyTime) >= 0.5]
  xTicks = xTicks[abs(xTicks - lastBiomarkerTime) >= 0.5]
  xTicks = c(xTicks, lastBiopsyTime, lastBiomarkerTime)
  xLabels = round(xTicks,2)
  xLabels[length(xLabels)-1] = paste0("t = ", round(lastBiopsyTime,2), " (Latest biopsy)")
  xLabels[length(xLabels)] = paste0("s = ", round(lastBiomarkerTime,2), " (Current visit)")
  xLabels = str_wrap(xLabels, width = 8)
  
  xLabelColors = rep(c("gray40", "black"), c(length(xLabels)-2,2))
  
  riskAxisBreaks = c(seq(0.25, 1,length.out = 4))
  riskAxisBreaks = c(0,riskAxisBreaks[abs(riskAxisBreaks - threshold) >= 0.1], threshold)
  riskAxisLabels = str_wrap(c(paste0(riskAxisBreaks[-length(riskAxisBreaks)] * 100, "%"), 
             paste0("k = ",threshold*100, "% (Biopsy threshold)")), width=8)
  riskAxisLabelColors = c(rep("gray40", length(riskAxisBreaks)-1), "black")
  
  if(usecolor == F){
    p=ggplot() +
      geom_line(data = psaDs, aes(x = visitTimeYears, y=fitted_log2psaplus1, linetype="Fitted PSA"), color="black") +
      geom_point(data = psaDs, size=POINT_SIZE, aes(x = visitTimeYears, y=log2psaplus1, shape="Observed PSA"), color="gray30") +
      geom_line(data = patientDs, aes(x = visitTimeYears, y=fitted_high_dre_prob, linetype="Fitted DRE"), color="black") +
      geom_point(data = patientDs, size=POINT_SIZE, aes(x = visitTimeYears, y=high_dre, shape="Observed DRE"), color="gray30") +
      geom_line(data = riskDf, aes(x=times, y=Mean), color="black") +
      geom_ribbon(data = riskDf, aes(x=times, ymin=Lower, ymax=Upper), 
                  fill="grey", alpha=0.5) +
      geom_label(aes(x=lastBiomarkerTime, y=maxMeanRisk, 
                     label=maxMeanRiskLabel, vjust = -0.25,hjust=0.75),
                     size = LABEL_SIZE) + 
      geom_vline(xintercept = lastBiopsyTime, linetype="solid") +
      #geom_vline(xintercept = lastBiomarkerTime, linetype="twodash") + 
      geom_segment(aes(x=-Inf, xend=lastBiopsyTime, y=maxYleft/2, yend=maxYleft/2), 
                   linetype="solid", color="gray", size=1) + 
      xlab("Follow-up time (years)") + 
      ylab(expression('Pr (DRE > T1c)            '*'log'[2]*'(PSA + 1)')) +
      scale_linetype_manual(name="",
                         labels= c("Fitted Pr (DRE > T1c)", expression('Fitted log'[2]*'(PSA + 1)')),
                         values = c("dotted", "dashed")) +       
      scale_shape_manual(name="",
                         labels=c("Observed DRE", expression('Observed log'[2]*'(PSA + 1)')),
                         values = c(17,16)) + 
      theme_bw() + 
      theme(text = element_text(size=FONT_SIZE), axis.text=element_text(size=FONT_SIZE),
            axis.line = element_line(),
            panel.grid.minor = element_blank(),
            axis.text.y = element_text(size=FONT_SIZE, color = rep(c("gray30", "gray30"), each=4)),
            axis.title.y = element_text(size=FONT_SIZE, color = "black"),
            axis.title.y.right = element_text(size=FONT_SIZE, color = "black"),
            axis.text.y.right = element_text(size=FONT_SIZE, color = riskAxisLabelColors),
            axis.text.x = element_text(size=FONT_SIZE, color = xLabelColors),
            legend.background = element_blank(), legend.position = "top",
            legend.text = element_text(size=FONT_SIZE-3))  +
      scale_x_continuous(breaks=xTicks, labels = xLabels) +
      scale_y_continuous(limits = c(minYLeft, maxYleft), 
                         breaks = c(seq(0, maxYleft/2 - DRE_PSA_Y_GAP, length.out = 4), seq(maxYleft/2 + DRE_PSA_Y_GAP, maxYleft, length.out = 4)),
                         labels = c(paste0(round(seq(0, 1, length.out = 4),2) * 100, "%"), 
                                    round(seq(0, maxPSA, length.out = 4),2)), 
                         sec.axis = sec_axis(~(.-minYLeft)/(maxYleft-minYLeft), 
                                             breaks= riskAxisBreaks,
                                             labels = riskAxisLabels,
                                             name = "Risk of cancer progression"))
  }else{
    p=NULL
  }
    
  return(p)
}

plotDynamicRiskProbNow = function(pid, fittedJointModel, maxVisitTime, 
                               maxPredictionTime=NA, lastBiopsyTime=NA, 
                               usecolor=F, threshold=0.15,
                               FONT_SIZE=12, POINT_SIZE = 2, DRE_PSA_Y_GAP=0.1,
                               LABEL_SIZE = 1){
  dataset = fittedJointModel$model_info$mvglmer_components$data
  
  patientDs = dataset[dataset$P_ID == pid & dataset$visitTimeYears<=maxVisitTime,]
  
  lastBiomarkerTime = max(patientDs$visitTimeYears)
  if(is.na(lastBiopsyTime)){
    lastBiopsyTime = getLastBiopsyTime(pid, upperLimitTime = maxVisitTime)
  }
  if(is.na(maxPredictionTime)){
    maxPredictionTime = lastBiomarkerTime
  }
  
  if(maxPredictionTime == lastBiopsyTime){
    stop("Max prediction time and last biopsy time are same")
  }
  
  futureTimes = seq(lastBiopsyTime, maxPredictionTime, length.out = 20)
  
  sfit = survfitJM(fittedJointModel, patientDs, idVar="P_ID", 
                   survTimes = futureTimes, last.time = lastBiopsyTime)
  
  patientDs$fitted_high_dre_prob = plogis(sfit$fitted.y[[1]]$high_dre)
  patientDs$fitted_log2psaplus1 = sfit$fitted.y[[1]]$log2psaplus1
  
  #The base of axes in this plot is of DRE
  minYLeft = 0
  maxYleft = 2 + DRE_PSA_Y_GAP * 2
  
  psaDs = patientDs[, c("visitTimeYears", "log2psaplus1", "fitted_log2psaplus1")]
  maxPSA = max(psaDs[,-1], na.rm=T)
  psaDs[,-1] = psaDs[,-1] / maxPSA + maxYleft/2 + DRE_PSA_Y_GAP
  
  riskDf = data.frame(sfit$summaries[[1]])
  riskDf = rbind(c("times"=futureTimes[1], "Mean"=1, "Median"=1, "Lower"=1, "Upper"=1), riskDf)
  # -1 because 1st column is time
  riskDf[,-1] = 1 - riskDf[,-1]
  maxMeanRiskLabel = paste0(round(max(riskDf$Mean)*100,2), "%")
  riskDf[,-1] = riskDf[, -1] * (maxYleft - minYLeft) + minYLeft
  maxMeanRisk = max(riskDf$Mean)
  
  xTicks = seq(0, maxPredictionTime, length.out = 6)
  xTicks = xTicks[abs(xTicks - lastBiopsyTime) >= 0.5]
  xTicks = xTicks[abs(xTicks - lastBiomarkerTime) >= 0.5]
  xTicks = c(xTicks, lastBiopsyTime, lastBiomarkerTime)
  xLabels = round(xTicks,2)
  xLabels[length(xLabels)-1] = paste0("t = ", round(lastBiopsyTime,2), " (Latest biopsy)")
  xLabels[length(xLabels)] = paste0("s = ", round(lastBiomarkerTime,2), " (Current visit)")
  xLabels = str_wrap(xLabels, width = 8)
  
  xLabelColors = rep(c("gray40", "black"), c(length(xLabels)-2,2))
  
  riskAxisBreaks = c(seq(0.25, 1,length.out = 4))
  riskAxisBreaks = c(0,riskAxisBreaks[abs(riskAxisBreaks - threshold) >= 0.1], threshold)
  riskAxisLabels = str_wrap(c(paste0(riskAxisBreaks[-length(riskAxisBreaks)] * 100, "%"), 
                              paste0("k = ",threshold*100, "% (Biopsy threshold)")), width=8)
  riskAxisLabelColors = c(rep("gray40", length(riskAxisBreaks)-1), "black")
  
  if(usecolor == F){
    p=ggplot() +
      geom_line(data = psaDs, aes(x = visitTimeYears, y=fitted_log2psaplus1, linetype="Fitted PSA"), color="black") +
      geom_point(data = psaDs, size=POINT_SIZE, aes(x = visitTimeYears, y=log2psaplus1, shape="Observed PSA"), color="gray30") +
      geom_line(data = patientDs, aes(x = visitTimeYears, y=fitted_high_dre_prob, linetype="Fitted DRE"), color="black") +
      geom_point(data = patientDs, size=POINT_SIZE, aes(x = visitTimeYears, y=high_dre, shape="Observed DRE"), color="gray30") +
      geom_line(data = riskDf, aes(x=times, y=Mean), color="black") +
      geom_ribbon(data = riskDf, aes(x=times, ymin=Lower, ymax=Upper), 
                  fill="grey", alpha=0.5) +
      geom_label(aes(x=lastBiomarkerTime, y=maxMeanRisk, 
                     label=maxMeanRiskLabel, vjust = -0.25,hjust=0.75),
                 size = LABEL_SIZE) + 
      geom_vline(xintercept = lastBiopsyTime, linetype="solid") +
      #geom_vline(xintercept = lastBiomarkerTime, linetype="twodash") + 
      geom_segment(aes(x=-Inf, xend=lastBiopsyTime, y=maxYleft/2, yend=maxYleft/2), 
                   linetype="solid", color="gray", size=1) + 
      xlab("Follow-up time (years)") + 
      ylab(expression('Pr (DRE > T1c)            '*'log'[2]*'(PSA + 1)')) +
      scale_linetype_manual(name="",
                            labels= c("Fitted Pr (DRE > T1c)", expression('Fitted log'[2]*'(PSA + 1)')),
                            values = c("dotted", "dashed")) +       
      scale_shape_manual(name="",
                         labels=c("Observed DRE", expression('Observed log'[2]*'(PSA + 1)')),
                         values = c(17,16)) + 
      theme_bw() + 
      theme(text = element_text(size=FONT_SIZE), axis.text=element_text(size=FONT_SIZE),
            axis.line = element_line(),
            panel.grid.minor = element_blank(),
            axis.text.y = element_text(size=FONT_SIZE, color = rep(c("gray30", "gray30"), each=4)),
            axis.title.y = element_text(size=FONT_SIZE, color = "black"),
            axis.title.y.right = element_text(size=FONT_SIZE, color = "black"),
            axis.text.y.right = element_text(size=FONT_SIZE, color = riskAxisLabelColors),
            axis.text.x = element_text(size=FONT_SIZE, color = xLabelColors),
            legend.background = element_blank(), legend.position = "top",
            legend.text = element_text(size=FONT_SIZE-3))  +
      scale_x_continuous(breaks=xTicks, labels = xLabels) +
      scale_y_continuous(limits = c(minYLeft, maxYleft), 
                         breaks = c(seq(0, maxYleft/2 - DRE_PSA_Y_GAP, length.out = 4), seq(maxYleft/2 + DRE_PSA_Y_GAP, maxYleft, length.out = 4)),
                         labels = c(paste0(round(seq(0, 1, length.out = 4),2) * 100, "%"), 
                                    round(seq(0, maxPSA, length.out = 4),2)), 
                         sec.axis = sec_axis(~(.-minYLeft)/(maxYleft-minYLeft), 
                                             breaks= riskAxisBreaks,
                                             labels = riskAxisLabels,
                                             name = "Risk of cancer progression"))
  }else{
    p=NULL
  }
  
  return(p)
}


plotObservedData = function(pid, fittedJointModel, maxVisitTime, 
                               lastBiopsyTime=NA, usecolor=F, 
                            FONT_SIZE=12, POINT_SIZE = 2, DRE_PSA_Y_GAP=0.1){
  dataset = fittedJointModel$model_info$mvglmer_components$data
  
  patientDs = dataset[dataset$P_ID == pid & dataset$visitTimeYears<=maxVisitTime,]
  
  lastBiomarkerTime = max(patientDs$visitTimeYears)
  if(is.na(lastBiopsyTime)){
    lastBiopsyTime = getLastBiopsyTime(pid, upperLimitTime = maxVisitTime)
  }
  
  secondLastBiopsyTime = getLastBiopsyTime(pid, upperLimitTime = maxVisitTime, lastnumber = 2)
  
  if(maxVisitTime == lastBiopsyTime){
    stop("Max prediction time and last biopsy time are same")
  }
  
  #The base of axes in this plot is of DRE
  minYLeft = 0
  maxYleft = 2 + DRE_PSA_Y_GAP * 2
  
  dreDs = patientDs[!is.na(patientDs$high_dre), c("visitTimeYears", "high_dre")]
  
  psaDs = patientDs[!is.na(patientDs$psa), c("visitTimeYears", "psa")]
  maxPSA = max(psaDs[,-1], na.rm=T)
  psaDs[,-1] = psaDs[,-1] / maxPSA + maxYleft/2 + DRE_PSA_Y_GAP
  
  fakeProgressionTime = ceiling(lastBiomarkerTime + 1)
  
  xTicks = seq(0, maxVisitTime, length.out = 6)
  xTicks = xTicks[abs(xTicks - lastBiopsyTime) >= 0.5]
  xTicks = xTicks[abs(xTicks - lastBiomarkerTime) >= 0.5]
  xTicks = xTicks[abs(xTicks - secondLastBiopsyTime) >= 0.5]
  xTicks = c(xTicks, fakeProgressionTime, secondLastBiopsyTime, lastBiopsyTime, lastBiomarkerTime)
  xLabels = round(xTicks,2)
  xLabels[length(xLabels)-3] = "\u221E"
  xLabels[length(xLabels)-2] = paste0(round(secondLastBiopsyTime,2), "\n(Older\nbiopsy)")
  xLabels[length(xLabels)-1] = paste0("t = ", round(lastBiopsyTime,2), "\n(Latest\nbiopsy)")
  xLabels[length(xLabels)] = paste0("s = ", round(lastBiomarkerTime,2), "\n(Current\nvisit)")
  
  xLabelColors = rep(c("gray40", "black"), c(length(xLabels)-4,4))
  
  if(usecolor == F){
    p=ggplot() +
      geom_line(data = psaDs, aes(x = visitTimeYears, y=psa), linetype="dashed") +
      geom_point(data = psaDs, size=POINT_SIZE, aes(x = visitTimeYears, y=psa, shape="Observed PSA"), color="gray30") +
      geom_line(data = dreDs, aes(x = visitTimeYears, y=high_dre), linetype="dashed") +
      geom_point(data = dreDs, size=POINT_SIZE, aes(x = visitTimeYears, y=high_dre, shape="Observed DRE"), color="gray30") +
      geom_vline(xintercept = lastBiopsyTime, linetype="solid") +
      geom_vline(xintercept = lastBiomarkerTime, linetype="twodash") + 
      geom_ribbon(aes(x=seq(lastBiopsyTime, fakeProgressionTime, length.out = 10), ymin=-Inf, ymax=Inf), fill="grey", alpha=0.5) + 
      geom_segment(aes(x=-Inf, xend=lastBiopsyTime, y=maxYleft/2, yend=maxYleft/2), 
                   linetype="solid", color="gray", size=1) + 
      geom_label(aes(x=lastBiomarkerTime, y=maxYleft/2, label="Biopsy now?")) + 
      xlab("Follow-up time (years)") + 
      ylab("DRE (binary)              PSA (ng/mL)") +
      scale_shape_manual(name="",
                         labels=c("Observed DRE (binary)", "Observed PSA (ng/mL)"),
                         values = c(17,16)) + 
      theme_bw() + 
      theme(text = element_text(size=FONT_SIZE), axis.text=element_text(size=FONT_SIZE),
            axis.line.x = element_line(arrow = arrow(length = unit(0.2,"cm"))),
            axis.line.y = element_line(),
            panel.grid.minor = element_blank(),
            axis.text.y = element_text(size=FONT_SIZE, color = rep(c("gray30", "gray30"), each=4)),
            axis.title.y = element_text(size=FONT_SIZE, color = "black"),
            axis.text.x = element_text(size=FONT_SIZE, color = xLabelColors),
            legend.background = element_blank(), legend.position = "top",
            legend.text = element_text(size=FONT_SIZE-3))  +
      scale_x_continuous(breaks=xTicks, labels = xLabels, limits = c(0, max(xTicks))) +
      scale_y_continuous(limits = c(minYLeft, maxYleft), 
                         breaks = c(seq(0, maxYleft/2 - DRE_PSA_Y_GAP, length.out = 2), seq(maxYleft/2 + DRE_PSA_Y_GAP, maxYleft, length.out = 4)), 
                         labels = str_wrap(c("T1c", "above T1c", round(seq(0, maxPSA, length.out = 4),2)), width=5))
  }else{
    p=NULL
  }
  
  return(p)
}

plotJMExplanationPlot = function(pid, fittedJointModel, maxVisitTime, 
                                 usecolor=F, FONT_SIZE=12, POINT_SIZE = 2, 
                                 DRE_PSA_Y_GAP=0.1, HAZARD_GAP = 1){
  
  generateTruePSAProfile = function(visitTimeYears, randomEff_psa){
    betas_psa = fittedJointModel$statistics$postMeans$betas2
    
    fixedPSAFormula = ~ 1 +I(Age - 70) +  I((Age - 70)^2) + ns(visitTimeYears, knots=c(0.1, 0.7, 4), Boundary.knots=c(0, 5.42))
    randomPSAFormula = ~ 1 + ns(visitTimeYears, knots=c(0.1, 0.7, 4), Boundary.knots=c(0, 5.42))
    
    df = data.frame(Age, visitTimeYears)
    model.matrix(fixedPSAFormula, df) %*% betas_psa + model.matrix(randomPSAFormula, df) %*% as.numeric(randomEff_psa)
  }
  
  generateTruePSASlope = function(visitTimeYears, randomEff_psa_slope){
    betas_psa_time = fittedJointModel$statistics$postMeans$betas2[4:7]
    
    fixedPSASlopeFormula = ~ 0 + dns(visitTimeYears, knots=c(0.1, 0.7, 4), Boundary.knots=c(0, 5.42))
    randomPSASlopeFormula = ~ 0 + dns(visitTimeYears, knots=c(0.1, 0.7, 4), Boundary.knots=c(0, 5.42))
    
    df = data.frame(visitTimeYears)
    model.matrix(fixedPSASlopeFormula, df) %*% betas_psa_time + model.matrix(randomPSASlopeFormula, df) %*% as.numeric(randomEff_psa_slope)
  }
  
  generateTrueDRELogOdds = function(visitTimeYears, randomEff_dre){
    betas_dre = fittedJointModel$statistics$postMeans$betas1
    
    fixedDREFormula = ~ 1 + I(Age - 70) +  I((Age - 70)^2) + visitTimeYears
    randomDREFormula = ~ 1 + visitTimeYears
    
    df = data.frame(Age, visitTimeYears)
    
    model.matrix(fixedDREFormula, df) %*% betas_dre + model.matrix(randomDREFormula, df) %*% as.numeric(randomEff_dre)
  }
  
  hazardFunc = function(visitTimeYears) {
    alphas = fittedJointModel$statistics$postMeans$alphas
    
    gammas = fittedJointModel$statistics$postMeans$gammas
    survivalFormula = ~ 0 + I(Age - 70) + I((Age - 70)^2)
    wGamma = as.numeric(model.matrix(survivalFormula, data = data.frame(Age=Age)) %*% gammas)
    
    baselinehazard = exp(splineDesign(fittedJointModel$control$knots, visitTimeYears, 
                                      ord = fittedJointModel$control$ordSpline, outer.ok = T) %*% fittedJointModel$statistics$postMeans$Bs_gammas)
    
    truePSA = generateTruePSAProfile(visitTimeYears, b_subject[3:7])
    truePSASlope = generateTruePSASlope(visitTimeYears, b_subject[4:7])
    trueDRELogOdds = generateTrueDRELogOdds(visitTimeYears, b_subject[1:2])
    
    y_Alpha = cbind(trueDRELogOdds, truePSA, truePSASlope) %*% alphas
    
    exp(wGamma + y_Alpha)
  }
  
  transformRange = function(data, newMin, newMax){
    oldMin = min(data, na.rm = T)
    oldMax = max(data, na.rm = T)
    
    (data - oldMin) * (newMax - newMin)/(oldMax-oldMin) + newMin
  }
  
  dataset = fittedJointModel$model_info$mvglmer_components$data
  patientDs = dataset[dataset$P_ID == pid & dataset$visitTimeYears<=maxVisitTime,]
  
  b_subject = mvJoint_dre_psa_dre_value$statistics$postMeans$b[which(unique(dataset$P_ID) == pid),]
  Age = patientDs$Age[1]
  
  dre_fixed_eff = mvJoint_dre_psa_dre_value$statistics$postMeans$betas1
  psa_fixed_eff = mvJoint_dre_psa_dre_value$statistics$postMeans$betas2
  
  times = seq(0, maxVisitTime, length.out = 20)
  truePSA = generateTruePSAProfile(times, b_subject[3:7])
  truePSASlope = generateTruePSASlope(times, b_subject[4:7])
  trueDREProbHighDRE = plogis(generateTrueDRELogOdds(times, b_subject[1:2]))
  trueHazard = hazardFunc(times)
  
  newMinDRE = 0
  newMaxDRE = 1/3 - DRE_PSA_Y_GAP/2
  newMinPSAVal = 1/3 + DRE_PSA_Y_GAP/2
  newMaxPSAVal = 2/3
  newMinPSASlope = 2/3 + DRE_PSA_Y_GAP
  newMaxPSASlope = 1 + DRE_PSA_Y_GAP/2
  
  transformedPSAVal = transformRange(c(truePSA, patientDs$log2psaplus1), newMinPSAVal, newMaxPSAVal)
  transformedTruePSASlope = transformRange(truePSASlope, newMinPSASlope, newMaxPSASlope)
  #I add a single 1 and remove it to make sure the range is taken correctly in DRE
  transformedHighDRE = transformRange(c(1, trueDREProbHighDRE, patientDs$high_dre), newMinDRE, newMaxDRE)[-1]
  transformedHazard = transformRange(trueHazard, 0, 1)
  
  patientDs$transformedLog2PSAPlus1 = transformedPSAVal[-c(1:length(truePSA))]
  patientDs$transformedHighDRE = transformedHighDRE[-c(1:length(trueDREProbHighDRE))]
  
  ggplot() + 
    geom_point(data = patientDs, size=POINT_SIZE, aes(x = visitTimeYears, y=transformedHighDRE, shape="Observed DRE"), color="gray30") +
    geom_point(data = patientDs, size=POINT_SIZE, aes(x = visitTimeYears, y=transformedLog2PSAPlus1, shape="Observed PSA"), color="gray30") +
    geom_line(aes(x = times, y=transformedPSAVal[1:length(truePSA)], linetype="Fitted PSA")) + 
    geom_line(aes(x = times, y=transformedTruePSASlope), linetype="solid") + 
    geom_line(aes(x = times, y=transformedHighDRE[1:length(trueDREProbHighDRE)], linetype="Fitted DRE")) + 
    geom_line(aes(x = times + maxVisitTime + HAZARD_GAP, y=transformedHazard)) +
    geom_text(aes(x = 4, y= newMinPSASlope, label="A")) + 
    geom_text(aes(x = 4, y= newMinPSAVal, label="B")) + 
    geom_text(aes(x = 4, y= newMinDRE + 0.1, label="C")) + 
    geom_text(aes(x = 9, y= newMinDRE, label="D")) + 
    
    geom_segment(aes(x=maxVisitTime+HAZARD_GAP/2, xend=maxVisitTime+HAZARD_GAP/2, y=-Inf, yend=Inf), 
                 linetype="solid", color="gray", size=1) + 
    geom_segment(aes(x=-Inf, xend=maxVisitTime+HAZARD_GAP/2, y=1/3, yend=1/3), 
                 linetype="solid", color="gray", size=1) +
    geom_segment(aes(x=-Inf, xend=maxVisitTime+HAZARD_GAP/2, y=2/3 + DRE_PSA_Y_GAP/2, yend=2/3+ DRE_PSA_Y_GAP/2), 
                 linetype="solid", color="gray", size=1) + 
    xlab("Follow-up time (years)                     Follow-up time (years)") + 
    ylab(expression('Pr (DRE > T1c)      '*'log'[2]*'(PSA + 1)'*'      Velocity: '*'log'[2]*'(PSA + 1)')) +
    scale_shape_manual(name="",
                       labels=c("Observed DRE (binary)", "Observed PSA (ng/mL)"),
                       values = c(17,16)) + 
    scale_linetype_manual(name="",
                          labels= c("Fitted Pr (DRE > T1c)", expression('Fitted log'[2]*'(PSA + 1)')),
                          values = c("dotted", "dashed")) +       
    theme_bw() + 
    theme(text = element_text(size=FONT_SIZE), axis.text=element_text(size=FONT_SIZE),
          axis.line = element_line(),
          panel.grid.minor = element_blank(),
          axis.text.y = element_text(size=FONT_SIZE, color = "gray30"),
          axis.title.y = element_text(size=FONT_SIZE, color = "black"),
          axis.text.y.right = element_text(size=FONT_SIZE, color = "black"),
          axis.text.x = element_text(size=FONT_SIZE, color=rep(c("gray40", "black"), each=4)),
          legend.background = element_blank(), legend.position = "top",
          legend.text = element_text(size=FONT_SIZE-3))  +
     scale_y_continuous(limits = c(0, newMaxPSASlope), 
                        breaks = c(seq(newMinDRE, newMaxDRE, length.out = 3),
                                   seq(newMinPSAVal, newMaxPSAVal, length.out = 3),
                                   seq(newMinPSASlope, newMaxPSASlope, length.out = 3)),
                        labels = c(c("0%", "50%", "100%"),
                                   round(seq(min(c(truePSA, patientDs$log2psaplus1), na.rm = T), max(c(truePSA, patientDs$log2psaplus1), na.rm = T), length.out = 3),2),
                                   round(seq(min(truePSASlope), max(truePSASlope), length.out = 3),2)),
                         sec.axis = sec_axis(~((.-min(trueHazard))*diff(range(trueHazard))/newMaxPSASlope + min(trueHazard)), 
                                             breaks= seq(min(trueHazard), max(trueHazard), length.out = 5),
                                             labels = round(seq(min(trueHazard), max(trueHazard), length.out = 5),2),
                                             name = "Hazard of cancer progression"))+
    scale_x_continuous(breaks = c(seq(0,maxVisitTime, length.out = 4),
                                  seq(maxVisitTime + HAZARD_GAP, maxVisitTime*2 + HAZARD_GAP, length.out = 4)),
                       labels = round(c(seq(0,maxVisitTime, length.out = 4), seq(0,maxVisitTime, length.out = 4)),2)) 
}

# obsDataPlot = plotObservedData(2340, mvJoint_dre_psa_dre_value, 4, FONT_SIZE = 11, DRE_PSA_Y_GAP = 0.25)
# ggsave(filename = "report/decision_analytic/mdm/latex/images/obsDataPlot_2340.eps",
#        plot=obsDataPlot, device=cairo_ps, height=4.5, width=6.1, dpi = 500)
# 
# dynRiskPlot = plotDynamicRiskProb(2340, mvJoint_dre_psa_dre_value, 4, FONT_SIZE = 11, DRE_PSA_Y_GAP = 0.25, LABEL_SIZE = 3)
# ggsave(filename = "report/decision_analytic/mdm/latex/images/obsDataPlot_2340_max4years.eps",
#        plot=dynRiskPlot, device=cairo_ps, height=4.5, width=6.1, dpi = 500)
# 
# ggsave(ggpubr::ggarrange(plotDynamicRiskProb(2340, mvJoint_dre_psa_dre_value, 4, FONT_SIZE = 11, DRE_PSA_Y_GAP = 0.25, LABEL_SIZE = 3), plotDynamicRiskProb(2340, mvJoint_dre_psa_dre_value, 6, FONT_SIZE = 11, DRE_PSA_Y_GAP = 0.25, LABEL_SIZE = 3), ncol = 1, nrow=2, common.legend = T), filename = "report/decision_analytic/mdm/latex/images/dynRiskPlot_2340.eps", width=6.1, height=7.5, device = cairo_ps)
# 
# jmExplanationPlot = plotJMExplanationPlot(1757, mvJoint_dre_psa_dre_value, 4, FONT_SIZE = 11)
# ggsave(filename = "report/decision_analytic/mdm/latex/images/jmExplanationPlot_1757.eps",
#        plot=jmExplanationPlot, device=cairo_ps, height=4.5, width=6.1, dpi = 500)
