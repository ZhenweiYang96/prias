load("Rdata/decision_analytic/DRE_PSA/mvJoint_dre_psa_dre_value.Rdata")
load("Rdata/decision_analytic/cleandata.Rdata")
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
                                      LABEL_SIZE = 1, xlim=NA){
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
                  fill="grey", alpha=0.4) +
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
                                  meanRiskProb = NA,
                                  threshold=0.10,
                                  FONT_SIZE=12, POINT_SIZE = 2, DRE_PSA_Y_GAP=0.1,
                                  LABEL_SIZE = 1, specialXticks = NA, xUpperlim = NA){
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
  
  sfit = survfitJM(fittedJointModel, patientDs, idVar="P_ID", 
                   survTimes = maxPredictionTime, last.time = lastBiopsyTime)
  
  patientDs$fitted_high_dre_prob = plogis(sfit$fitted.y[[1]]$high_dre)
  patientDs$fitted_log2psaplus1 = sfit$fitted.y[[1]]$log2psaplus1
  
  #The base of axes in this plot is of DRE
  minYLeft = 0
  maxYleft = 2 + DRE_PSA_Y_GAP * 2
  
  psaDs = patientDs[, c("visitTimeYears", "log2psaplus1", "fitted_log2psaplus1")]
  maxPSA = max(psaDs[,-1], na.rm=T)
  minPSA = min(psaDs[,-1], na.rm = T)
  
  newMinPSA = (maxYleft/2 + DRE_PSA_Y_GAP)
  newMaxPSA = maxYleft
  
  psaDs[,-1] = ((psaDs[,-1] - minPSA) * (newMaxPSA - newMinPSA))/ (maxPSA - minPSA) + newMinPSA
  
  if(is.na(meanRiskProb)){
    meanRiskProb = 1 - sfit$summaries[[1]][, "Mean"]
  }
  maxMeanRiskScaled = meanRiskProb * (maxYleft - minYLeft) + minYLeft
  thresholdScaled = threshold *  (maxYleft - minYLeft) + minYLeft
  
  
  if(is.na(xUpperlim)){
    xUpperlim = 1.1 * maxPredictionTime
  }
  
  col_width = 0.1 * xUpperlim
  
  if(is.na(specialXticks)){
    specialXticks = 0
  }else{
    specialXticks = c(0, specialXticks)
  }
  
  xTicks = c(specialXticks, lastBiopsyTime, maxPredictionTime)
  xLabels = c(specialXticks, paste0("t = ", round(lastBiopsyTime,1), "\n(Latest\nbiopsy)"),
              paste0("s = ", round(maxPredictionTime,1), "\n(Current\nvisit)"))
  
  riskFillColor = if(meanRiskProb < threshold){
    "forestgreen"
  }else{
    "red3"
  }
  
  p=ggplot() +
    geom_segment(aes(x=-Inf, xend=maxPredictionTime, y=maxYleft/2, yend=maxYleft/2), 
                 linetype="solid") + 
    geom_col(aes(x=c(maxPredictionTime,maxPredictionTime), y=c(maxMeanRiskScaled, maxYleft-maxMeanRiskScaled)), 
             fill=c(riskFillColor, "white"),  color=riskFillColor, width = col_width) + 
    geom_vline(xintercept = lastBiopsyTime, linetype="solid") +
    geom_line(data = psaDs, aes(x = visitTimeYears, y=fitted_log2psaplus1, linetype="Fitted PSA")) +
    geom_point(data = psaDs, size=POINT_SIZE, aes(x = visitTimeYears, y=log2psaplus1, shape="Observed PSA")) +
    geom_line(data = patientDs, aes(x = visitTimeYears, y=fitted_high_dre_prob, linetype="Fitted DRE")) +
    geom_point(data = patientDs, size=POINT_SIZE, aes(x = visitTimeYears, y=high_dre, shape="Observed DRE")) +
    geom_segment(aes(x=maxPredictionTime-col_width/2, xend=maxPredictionTime + col_width, y=thresholdScaled, yend=thresholdScaled), 
                 linetype="dashed", color="firebrick1", size=0.5) +
    geom_label(aes(x=maxPredictionTime, y=1.6, 
                   label=paste0("Cancer progression\nrisk at current visit:\n", round(meanRiskProb * 100,2), "%")),
               color = riskFillColor, size=LABEL_SIZE)+
    geom_text(aes(x=maxPredictionTime + col_width*1.5, y=thresholdScaled, 
                  label=paste0("k = ",threshold*100, "%\n(Biopsy threshold)")),
              color = "firebrick1", size=LABEL_SIZE)+
    xlab("Follow-up time (years)") + 
    ylab(expression('Pr (DRE > T1c)            '*'log'[2]*'(PSA + 1)')) +
    scale_linetype_manual(name="",
                          labels= c("Fitted Pr (DRE > T1c)", expression('Fitted log'[2]*'(PSA + 1)')),
                          values = c("dotted", "dashed")) +       
    scale_shape_manual(name="",
                       labels=c("Observed DRE\n(T1c / above T1c)", expression('Observed log'[2]*'(PSA + 1)')),
                       values = c(17,16)) + 
    theme_bw() + 
    theme(text = element_text(size=FONT_SIZE), axis.text=element_text(size=FONT_SIZE),
          axis.line = element_line(),
          
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          legend.background = element_blank(), legend.position = "bottom",
          legend.direction = "horizontal",
          legend.text = element_text(size=FONT_SIZE-3)) +
    scale_x_continuous(breaks=xTicks, labels = xLabels, limits = c(0, xUpperlim)) +
    scale_y_continuous(limits = c(minYLeft, maxYleft), 
                       breaks = c(seq(0, maxYleft/2 - DRE_PSA_Y_GAP, length.out = 3), seq(maxYleft/2 + DRE_PSA_Y_GAP, maxYleft, length.out = 5)),
                       labels = c(paste0(round(seq(0, 1, length.out = 3),2) * 100, "%"), 
                                  round(seq(minPSA, maxPSA, length.out = 5),2)))
  
  return(p)
}

plotObservedData = function(pid, fittedJointModel, maxVisitTime, 
                            lastBiopsyTime=NA, 
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
  
  xTicks = c(0, lastBiopsyTime, lastBiomarkerTime, fakeProgressionTime)
  xLabels = c(0, paste0("t = ", round(lastBiopsyTime,1), "\n(Latest\nbiopsy)"),
              paste0("s = ", round(lastBiomarkerTime,1), "\n(Current\nvisit)"),
              "\u221E")
  
  p=ggplot() +
    geom_point(data = psaDs, size=POINT_SIZE, aes(x = visitTimeYears, y=psa, shape="Observed PSA")) +
    geom_point(data = dreDs, size=POINT_SIZE, aes(x = visitTimeYears, y=high_dre, shape="Observed DRE")) +
    geom_vline(xintercept = lastBiopsyTime, linetype="solid") +
    geom_vline(xintercept = lastBiomarkerTime, linetype="twodash") + 
    geom_ribbon(aes(x=seq(lastBiopsyTime, fakeProgressionTime, length.out = 10), 
                    ymin=-Inf, ymax=Inf, fill="Region with risk of\n cancer progression"), alpha=0.2) + 
    geom_segment(aes(x=-Inf, xend=lastBiomarkerTime, y=maxYleft/2, yend=maxYleft/2), 
                 linetype="solid") + 
    geom_label(aes(x=lastBiomarkerTime, y=maxYleft/2, label="Biopsy\nnow?"), color="firebrick1") + 
    xlab("Follow-up time (years)") + 
    ylab("DRE (binary)               PSA (ng/mL)") +
    scale_fill_manual(name="", values="red2")+
    scale_shape_manual(name="",
                       labels=c("Observed DRE\n(T1c / above T1c)", "Observed PSA\n(ng/mL)"),
                       values = c(17,16)) + 
    theme_bw() + 
    theme(text = element_text(size=FONT_SIZE), axis.text=element_text(size=FONT_SIZE),
          axis.line.x = element_line(arrow = arrow(length = unit(0.2,"cm"))),
          axis.line.y = element_line(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          legend.background = element_blank(), legend.position = "bottom",
          legend.text = element_text(size=FONT_SIZE-3),
          plot.margin = margin(0, 0, 0, 0, "pt")) +
    scale_x_continuous(breaks=xTicks, labels = xLabels, limits = c(0, max(xTicks))) +
    scale_y_continuous(limits = c(minYLeft, maxYleft), 
                       breaks = c(seq(0, maxYleft/2 - DRE_PSA_Y_GAP, length.out = 2), seq(maxYleft/2 + DRE_PSA_Y_GAP, maxYleft, length.out = 5)), 
                       labels = c("T1c", "above\nT1c", round(seq(0, maxPSA, length.out = 5),1)))
  
  return(p)
}

plotJMExplanationPlot_Stacked = function(pid, fittedJointModel, maxVisitTime, 
                                         usecolor=F, FONT_SIZE=12, POINT_SIZE = 2){
  
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
  
  dataset = fittedJointModel$model_info$mvglmer_components$data
  patientDs = dataset[dataset$P_ID == pid & dataset$visitTimeYears<=maxVisitTime,]
  
  b_subject = fittedJointModel$statistics$postMeans$b[which(unique(dataset$P_ID) == pid),]
  Age = patientDs$Age[1]
  
  dre_fixed_eff = fittedJointModel$statistics$postMeans$betas1
  psa_fixed_eff = fittedJointModel$statistics$postMeans$betas2
  
  times = seq(0, maxVisitTime, length.out = 20)
  truePSA = generateTruePSAProfile(times, b_subject[3:7])
  truePSASlope = generateTruePSASlope(times, b_subject[4:7])
  trueDREProbHighDRE = plogis(generateTrueDRELogOdds(times, b_subject[1:2]))
  trueHazard = hazardFunc(times)
  
  hazardPlotYbreaks = seq(min(trueHazard), 
                          max(trueHazard), length.out = 3)
  hazardPlot = ggplot() + geom_line(aes(x=times, y=trueHazard), color='red3') + 
    theme_bw() + 
    ylab("Estimated Hazard of\ncancer progression") + xlab("Follow-up time (years)") +
    theme(text = element_text(size=FONT_SIZE), 
          axis.title.y = element_text(colour = 'red3'),
          axis.text.y = element_text(colour = 'red3', size=FONT_SIZE),
          axis.line = element_line(),
          axis.text=element_text(size=FONT_SIZE),
          plot.margin = margin(0, 5.5, 0, 5.5, "pt")) +
    scale_y_continuous(breaks = hazardPlotYbreaks, labels = round(hazardPlotYbreaks,1))
  
  
  #also added fake points and lines in the DRE plot to add legend for points in ggpubr
  drePlot = ggplot() + 
    geom_point(data = patientDs, size=POINT_SIZE, 
               aes(x = visitTimeYears, y=high_dre, shape="Observed DRE")) +
    geom_point(aes(x = 0, y=10, shape="Observed PSA")) +
    geom_line(aes(x = times, y=trueDREProbHighDRE, linetype="Fitted DRE")) +
    geom_line(aes(x = 0, y=10, linetype="Fitted PSA")) +
    scale_shape_manual(name="",
                       labels=c("Observed DRE\n(T1c / above T1c)", expression('Observed log'[2]*'(PSA + 1)')),
                       values = c(17, 16)) + 
    scale_linetype_manual(name="",
                          labels= c("Fitted Pr (DRE > T1c)", expression('Fitted log'[2]*'(PSA + 1)')),
                          values = c("dotted", "dashed")) +       
    theme_bw() + ylab("Probability (DRE > T1c)") + 
    scale_y_continuous(breaks = seq(0,1, length.out = 3), 
                       labels=paste0(seq(0,1,length.out = 3)*100, "%"), limits = c(0,1),
                       sec.axis = dup_axis(trans = ~.,
                                           breaks = c(0,1),
                                           labels = c("T1c", "\nabove\nT1c"),
                                           name = "Observed DRE (binary)"))+
    theme(text = element_text(size=FONT_SIZE), 
          axis.line = element_line(),
          axis.ticks.x = element_blank(), 
          axis.title.x = element_blank(), 
          axis.text.x = element_blank(),
          legend.text = element_text(size = FONT_SIZE-3),
          plot.margin = margin(2, 5.5, 0, 5.5, "pt"))
  
  psaPlotYbreaks = seq(min(c(truePSA,patientDs$log2psaplus1), na.rm = T),
                       max(c(truePSA,patientDs$log2psaplus1), na.rm = T),
                       length.out = 3)
  
  psaPlot = ggplot() + 
    geom_point(data = patientDs, size=POINT_SIZE, 
               aes(x = visitTimeYears, y=log2psaplus1), 
               shape=16) +
    geom_line(aes(x = times, y=truePSA), linetype="dashed") +
    theme_bw() + ylab(expression('log'[2]*'(PSA + 1) Levels')) + 
    theme(text = element_text(size=FONT_SIZE), 
          axis.line = element_line(),
          axis.ticks.x = element_blank(), 
          axis.title.x = element_blank(), 
          axis.text.x = element_blank(), 
          plot.margin = margin(0, 5.5, 0, 5.5, "pt")) +
    scale_y_continuous(breaks = psaPlotYbreaks, labels = round(psaPlotYbreaks,1))
  
  psaVelocityPlotYbreaks = seq(min(truePSASlope), 
                               max(truePSASlope), length.out = 3)
  
  psaVelocityPlot = ggplot() + 
    geom_line(aes(x = times, y=truePSASlope)) +
    theme_bw() + ylab(bquote(atop("Estimated", 'log'[2]*'(PSA + 1) Velocity'))) + 
    theme(text = element_text(size=FONT_SIZE), 
          axis.ticks.x = element_blank(), 
          axis.title.x = element_blank(), 
          axis.text.x = element_blank(), 
          axis.line = element_line(),
          plot.margin = margin(0, 5.5, 0, 5.5, "pt"))+
    scale_y_continuous(breaks = psaVelocityPlotYbreaks, 
                       labels = round(psaVelocityPlotYbreaks,1))
  
  ggpubr::ggarrange(drePlot, psaPlot, psaVelocityPlot, hazardPlot, labels = "AUTO", 
                    hjust = -10, common.legend = T, legend = "bottom", heights = c(1,1,1,1.20),
                    nrow=4,ncol=1, align = "v")
}

threshold_choice_plot = function(f1threshold="7.5%",lastBiopsyTime=1, currentVisitTime=2,
                                 intermediateXMarker = 3,
                                 maxTime=4, FONT_SIZE=12){
  
  yfixed = 3
  yf1 = 2
  
  ggplot() + geom_ribbon(aes(x=c(lastBiopsyTime, maxTime),ymin=-Inf, ymax=Inf, 
                             fill="Region with risk of cancer progression"), alpha=0.2) +
    geom_vline(xintercept = lastBiopsyTime) + 
    geom_vline(xintercept = currentVisitTime, linetype="twodash") +
    geom_segment(aes(x=-Inf, xend=currentVisitTime, y=0.5*(yf1+yfixed), yend=0.5*(yf1+yfixed)))+
    geom_segment(aes(x=lastBiopsyTime, xend=currentVisitTime, y=0.25, yend=0.25), 
                 arrow = arrow(length = unit(0.2, "cm"), ends="both", type="closed")) +
    geom_label(aes(x=currentVisitTime, y=yfixed, label="Fixed\nthreshold = 10%"),color='red3') + 
    geom_label(aes(x=currentVisitTime, y=yf1, 
                   label=paste0("Dynamic\nthreshold = ", f1threshold)), color='red3') + 
    geom_text(aes(x=0.5*(currentVisitTime+lastBiopsyTime),y=3, 
                  label="Common threshold\n at all visits.")) + 
    geom_text(aes(x=0.5*(currentVisitTime+lastBiopsyTime),y=1, 
                  label="Threshold is chosen\nfor a combination of\nlatest biopsy time 't' and\ntime gap since last biopsy 's-t',\nby maximizing\n F1 score."))+
    scale_fill_manual(name="", values="red2")+
    scale_x_continuous(breaks=c(0,lastBiopsyTime,currentVisitTime,intermediateXMarker, maxTime), 
                       labels=c("0", paste0("t=", lastBiopsyTime, "\n(Latest\nBiopsy)"), 
                                paste0("s=",currentVisitTime, "\n(Current\nvisit)"), 
                                paste(intermediateXMarker), "\u221E"), limits = c(0,maxTime)) + 
    scale_y_continuous(breaks=c(yf1,yfixed), 
                       labels=c("Threshold based \non F1 score", 
                                "Fixed threshold\non all visits"),limits = c(0,yfixed+0.25)) + 
    theme_bw() + theme(text = element_text(size=FONT_SIZE), 
                       axis.line = element_line(), legend.position = "bottom")+
    xlab("Follow-up time (years)") + ylab("Choice of risk threshold") 
}

obsDataPlot = plotObservedData(2340, mvJoint_dre_psa_dre_value, 4, FONT_SIZE = 11, 
                                POINT_SIZE = 3, DRE_PSA_Y_GAP = 0.3)
ggsave(filename = "report/decision_analytic/mdm/latex/images/obsDataPlot_2340.eps",
        plot=obsDataPlot, device=cairo_ps, height=5.5/1.333, width=5.5, dpi = 500)
# # 
# # dynRiskPlot = plotDynamicRiskProb(2340, mvJoint_dre_psa_dre_value, 4, FONT_SIZE = 11, DRE_PSA_Y_GAP = 0.25, LABEL_SIZE = 3)
# # ggsave(filename = "report/decision_analytic/mdm/latex/images/obsDataPlot_2340_max4years.eps",
# #        plot=dynRiskPlot, device=cairo_ps, height=4.5, width=6.1, dpi = 500)
# # 
dynRiskPlot1 = plotDynamicRiskProbNow(2340, mvJoint_dre_psa_dre_value, 4, lastBiopsyTime = 2.56, 
                                      meanRiskProb = 0.078, FONT_SIZE = 11, DRE_PSA_Y_GAP = 0.25, LABEL_SIZE = 3, 
                                      xUpperlim = 6.5, POINT_SIZE = 3)
dynRiskPlot1 = dynRiskPlot1 + ggtitle("Biopsy not recommended for patient j at year 4") + 
   theme(plot.title = element_text(color = "forestgreen"))
#plot.margin = margin(t=5,b=30))
dynRiskPlot2 = plotDynamicRiskProbNow(2340, mvJoint_dre_psa_dre_value, 8, lastBiopsyTime = 2.56, meanRiskProb = 0.135, FONT_SIZE = 11, DRE_PSA_Y_GAP = 0.25, LABEL_SIZE = 3, specialXticks = 4,
                                      xUpperlim = 6.5, POINT_SIZE = 3) 
dynRiskPlot2 = dynRiskPlot2  + ggtitle("Biopsy recommended for patient j at year 5.3") + 
 theme(plot.title = element_text(color="red3"), legend.position = "none")
ggsave(ggpubr::ggarrange(dynRiskPlot1, dynRiskPlot2,
                        ncol = 1, nrow=2, labels = "AUTO", align = "v", heights = c(1.15,1)), filename = "report/decision_analytic/mdm/latex/images/dynRiskPlot_2340.eps", 
        width=7, height=9, device = cairo_ps, dpi=500)
#  
# # jmExplanationPlot = plotJMExplanationPlot(1757, mvJoint_dre_psa_dre_value, 4, FONT_SIZE = 11)
jmExplanationPlot = plotJMExplanationPlot_Stacked(113, mvJoint_dre_psa_dre_value, 4,
                                                   POINT_SIZE = 3, FONT_SIZE = 12)
ggsave(filename = "report/decision_analytic/mdm/latex/images/jmExplanationPlot_1757.eps",
        plot=jmExplanationPlot, device=cairo_ps, height=8.5, width=7, dpi = 500)

###Threshold choice plot
threshold_choice_plot1 = threshold_choice_plot(f1threshold = "8.7%", currentVisitTime = 3, intermediateXMarker = 3.5, maxTime = 4) + 
  ggtitle("Latest biopsy time t=1 year, and current visit time s=3 years")
threshold_choice_plot2 = threshold_choice_plot(f1threshold = "5.9%",currentVisitTime = 3.5, intermediateXMarker = NA, maxTime = 4) + 
  ggtitle("Latest biopsy time t=1 year, and current visit time s=3.5 years") + theme(legend.position = "none")
ggsave(ggpubr::ggarrange(threshold_choice_plot1, threshold_choice_plot2,
                         ncol = 1, nrow=2, labels = "AUTO", align = "v", heights = c(1.15,1)),
       filename = "report/decision_analytic/mdm/latex/images/threshold_choice_plot.eps", 
       width=7, height=9, device = cairo_ps, dpi=500)
