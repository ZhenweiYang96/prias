#load("~/Google Drive/PhD/src/prias/Rdata/Gleason as event/tdist/log2psa_plus1_andplus_pt1/mvjoint_log2psa_plus1_spline_pt1pt54_pt1.Rdata")
#training_psa_data_set=mvjoint_log2psa_plus1_spline_pt1pt54_pt1$model_info$mvglmer_components$data
#training.prias.id=mvjoint_log2psa_plus1_spline_pt1pt54_pt1$model_info$coxph_components$data
#load("~/Google Drive/PhD/src/prias/Rdata/Gleason as event/cleandata.Rdata")



calculateHazard = function(mvJointModelObject, patientId, maxRowNum=NA){
  
  fittedObjIndex = which(training.prias.id$P_ID==patientId)
  fittedRandEff = apply(mvJointModelObject$mcmc$b, c(1,2), mean)[fittedObjIndex, ]
  
  ds = training_psa_data_set[training_psa_data_set$P_ID == patientId, ] 
  Age = ds$Age[1]
  
  lastVisitTime = max(ds$visitTimeYears)
  lastBiopsyTime = tail(ds$visitTimeYears[!is.na(ds$gleason)],1)
  
  baselinehazard = exp(splineDesign(mvJointModelObject$control$knots, ds$visitTimeYears, 
                                ord = mvJointModelObject$control$ordSpline, 
                                outer.ok = T) %*% mvJointModelObject$statistics$postMeans$Bs_gammas)
  
  fixedValueFormula = ~ 1 + I(Age - 70) + I((Age - 70)^2) + ns(visitTimeYears, knots = c(0.1, 0.5, 4), Boundary.knots = c(0, 7))
  randomValueFormula = ~ 1 + ns(visitTimeYears, knots = c(0.1), Boundary.knots = c(0, 7))
  
  fixedSlopeFormula = ~ 0 + dns(visitTimeYears, knots = c(0.1, 0.5, 4), Boundary.knots = c(0, 7))
  randomSlopeFormula = ~ 0 + dns(visitTimeYears, knots = c(0.1), Boundary.knots = c(0, 7))
  
  survivalFormula = ~ 0 + I(Age - 70) + I((Age - 70)^2)
  
  fittedValue = model.matrix(fixedValueFormula, data=ds) %*% mvJointModelObject$statistics$postMeans$betas1 +
    model.matrix(randomValueFormula, data=ds) %*% fittedRandEff
    
  fittedVelocity = model.matrix(fixedSlopeFormula, data=ds) %*% mvJointModelObject$statistics$postMeans$betas1[4:7] +
    model.matrix(randomSlopeFormula, data=ds) %*% fittedRandEff[-1]
  
  fittedBaselineCov_in_hazard = model.matrix(survivalFormula, data=ds) %*% mvJointModelObject$statistics$postMeans$gammas
  
  nonBaseLinehazard = exp(fittedBaselineCov_in_hazard + 
                            mvJointModelObject$statistics$postMeans$alphas[1] * fittedValue +
                            mvJointModelObject$statistics$postMeans$alphas[2] * fittedVelocity)
  
  totalHazard = baselinehazard * nonBaseLinehazard
}


plotLog2PSAPlus1JMFit = function(mvJointModelObject, patientId, maxRowNum=NA){
  
  ds = prias_long[prias_long$P_ID == patientId, ] 
  ds$log2psa_plus1 = log(ds$psa + 1, 2)
  Age = ds$Age[1]
  
  if(!is.na(maxRowNum)){
    ds = ds[1:maxRowNum,]  
  }
  
  lastVisitTime = max(ds$visitTimeYears)
  lastBiopsyTime = tail(ds$visitTimeYears[!is.na(ds$gleason)],1)

  fittedDf = data.frame("Age"=Age,"visitTimeYears"=ds$visitTimeYears, 
                        "log2psa_plus1"=survfitJM(mvJointModelObject, ds, idVar="P_ID", last.time = lastBiopsyTime)$fitted.y[[1]][[1]], 
                        type="Fitted")
  
  ds = ds[,c("Age", "visitTimeYears", "log2psa_plus1")]
  ds$type="Observed"
  ds = rbind(ds, fittedDf)
  
  breaks = if(round(max(ds$visitTimeYears))>6){
    seq(0, round(max(ds$visitTimeYears)), by=2)
  }else if(round(max(ds$visitTimeYears))>2){
    seq(0, round(max(ds$visitTimeYears)), by=1)
  }else if(round(max(ds$visitTimeYears))==2){
    seq(0, round(max(ds$visitTimeYears)), length.out=5)
  }else{
    seq(0, round(max(ds$visitTimeYears)), length.out=3)
  }

  linetypes = c("dashed", "solid", "solid")  
  
  p=ggplot(data=ds[!is.na(ds$log2psa_plus1),]) + 
    geom_line(aes(x=visitTimeYears, y=log2psa_plus1, linetype=type)) + 
    geom_point(size=2, data=ds[ds$type=="Observed",], aes(x=visitTimeYears, y=log2psa_plus1)) +
    scale_x_continuous(breaks=breaks) + 
    geom_vline(aes(xintercept = lastBiopsyTime, linetype="Latest biopsy"), show.legend=F) + 
    scale_linetype_manual(values=linetypes) + 
    theme(text = element_text(size=11), axis.text=element_text(size=11),
          legend.text = element_text(size=11), legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5, size=13)) +
    xlab("Time (Years)") + ylab(expression('log'[2]*' (PSA + 1)')) + ggtitle(paste("Last visit time =", round(lastVisitTime, 2), "years"))
  
  print(p)
  
  return(p)
}

plotMarginalPSAFittedCurveJM = function(mvJointModelObject){
  newDF <- data.frame(visitTimeYears = seq(0, 10, by = 0.5), Age = 70)
  
  Xmatrix = model.matrix(~ 1 + I(Age - 70) + I((Age - 70)^2) + ns(visitTimeYears, knots = c(0.1, 0.5, 4), Boundary.knots = c(0, 7)), newDF)
  fitted = Xmatrix %*% t(mvjoint_log2psa_plus1_spline_pt1pt54_pt1$mcmc$betas1)
  newDF$fittedMean = apply(fitted, 1, mean)
  newDF$fittedLowCI = apply(fitted, 1, quantile, probs=0.025)
  newDF$fittedUpCI = apply(fitted, 1, quantile, probs=0.975)
  
  plot = ggplot(data=newDF) + geom_ribbon(aes(x=visitTimeYears, ymin=fittedLowCI, ymax=fittedUpCI), fill="grey80") + 
    geom_line(aes(y=fittedMean, x=visitTimeYears)) + 
    scale_x_continuous(breaks=seq(0, 10, by = 1)) + 
    scale_y_continuous(breaks=seq(0, 25, by = 0.1)) + 
    xlab("Time (Years)") + 
    ylab(expression('Fitted marginal '*'log'[2]*' (PSA + 1)')) + 
    theme(text = element_text(size=11), axis.text=element_text(size=11), 
          legend.position = "top", legend.title = element_blank(), 
          legend.text = element_text(size=11))
  
  return(plot)
}


getFittedVelocities = function(mvJointModelObject){ 
  
  df = mvJointModelObject$model_info$mvglmer_components$data
  df$P_ID = droplevels(df$P_ID)
  
  betas = mvJointModelObject$statistics$postMeans$betas1
  b = mvJointModelObject$statistics$postMeans$b
  
  fixedSlopeFormula = ~ 0 + dns(visitTimeYears, knots = c(0.1, 0.5, 4), Boundary.knots = c(0, 7))
  randomSlopeFormula = ~ 0 + dns(visitTimeYears, knots = c(0.1), Boundary.knots = c(0, 7))
  
  Xbeta = model.matrix(fixedSlopeFormula,data=df) %*% betas[-c(1:3)]
  Zb = rowSums(model.matrix(randomSlopeFormula, data=df) * b[as.numeric(df$P_ID),-1])
  
  return(Xbeta + Zb)
}

#fitted profiles for 911
# profilesfit911 = lapply(12:17, plotLog2PSAPlus1JMFit, mvJointModelObject=mvjoint_log2psa_plus1_spline_pt1pt54_pt1, patientId=911)
# ggsave(file="report/pers_sched/latex/biometrics_submission/images/model_fit/fitted_911_log2psaplus1.eps", 
#        grid_arrange_shared_legend(plotlist = profilesfit911, ncol=2, nrow=3),
#        width=8.27, height=8.27)
# 
# profilesfit2340 = lapply(12:17, plotLog2PSAPlus1JMFit, mvJointModelObject=mvjoint_log2psa_plus1_spline_pt1pt54_pt1, patientId=2340)
# ggsave(file="report/pers_sched/latex/biometrics_submission/images/model_fit/fitted_2340_log2psaplus1.eps", 
#        grid_arrange_shared_legend(plotlist = profilesfit2340, ncol=2, nrow=3),
#        width=8.27, height=8.27)
# 
# profilesfit3174 = lapply(4:9, plotLog2PSAPlus1JMFit, mvJointModelObject=mvjoint_log2psa_plus1_spline_pt1pt54_pt1, patientId=3174)
# ggsave(file="report/pers_sched/latex/biometrics_submission/images/model_fit/fitted_3174_log2psaplus1.eps", 
#        grid_arrange_shared_legend(plotlist = profilesfit3174, ncol=2, nrow=3),
#        width=8.27, height=8.27)
