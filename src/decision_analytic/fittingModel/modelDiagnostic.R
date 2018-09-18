source("src/decision_analytic/fittingModel/weighted_fitted.R")

ggsave(plotFittedDREMarginal(mvJoint_dre_psa_dre_value, FONT_SIZE = 9), file="report/decision_analytic/mdm/latex/images/marginal_dre.eps", device = cairo_ps, dpi = 500, height=4.5/1.3333, width=4.5)
#plots marginal prob by default, otherwise if logOdds = T, then plots marginal log odds
plotFittedDREMarginal = function(modelObject, logOdds = F, FONT_SIZE=10){
  visitTimeYears = seq(0, 10, 0.5)
  df = data.frame(Age = 70, visitTimeYears)
  
  betas_dre_matrix = modelObject$mcmc$betas1
  D_matrix = modelObject$mcmc$D
  
  fixedDREFormula = ~ 1 + I(Age - 70) +  I((Age - 70)^2) + visitTimeYears
  randomDREFormula = ~ 1 + visitTimeYears
  
  totalMCMC = dim(D_matrix)[1]
  
  xbetaZb_mcmc = sapply(1:totalMCMC, function(m){
    beta_m = betas_dre_matrix[m,]
    D_m = D_matrix[m,,]
    
    nSub = 10000
    b_m = MASS::mvrnorm(nSub, mu = rep(0, nrow(D_m)), D_m)[,1:2]
    
    xbeta = model.matrix(fixedDREFormula, df) %*% beta_m
    Zb = model.matrix(randomDREFormula, df) %*% t(b_m)
    
    xbetaZb = matrix(rep(c(xbeta), nSub), ncol=nSub, byrow = F) + Zb
    apply(plogis(xbetaZb), 1, mean)
  })
  
  if(logOdds==T){
    xbetaZb_mcmc = log(xbetaZb_mcmc / (1-xbetaZb_mcmc))
  }
  
  meanXbetaZb = apply(xbetaZb_mcmc, 1, mean)
  upperXbetaZb = apply(xbetaZb_mcmc, 1, quantile, probs=0.975)
  lowerXbetaZb = apply(xbetaZb_mcmc, 1, quantile, probs=0.025)
  
  plot = ggplot() + 
    geom_line(aes(x=visitTimeYears, y=meanXbetaZb)) +
    geom_ribbon(aes(x=visitTimeYears, ymin = lowerXbetaZb, ymax = upperXbetaZb), fill="grey", alpha=0.5) + 
    theme_bw() + 
    theme(text = element_text(size=FONT_SIZE), axis.text=element_text(size=FONT_SIZE),
          axis.line = element_line(),
          panel.grid.minor = element_blank()) + 
    xlab("Follow-up time (years)")
  
  if(logOdds==F){
   plot = plot + scale_y_continuous(breaks = seq(0,1, by = 0.25), labels=paste0(seq(0,1,by=0.25)*100, "%"),limits = c(0,1)) +
    ylab(expression('Probability (DRE > T1c)'))
  }else{
    plot = plot + ylab(expression('log odds (DRE > T1c)'))
  }
}

#set.seed(2018)
#temp = lapply(1:9, plotFittedDRESubject, modelObject=mvJoint_dre_psa_dre_value, FONT_SIZE=13, showTitle=F)
#subjectplot = ggpubr::ggarrange(plotlist = temp, nrow = 3, ncol=3, common.legend = T, legend = "bottom")
#ggsave(subjectplot, file="report/decision_analytic/mdm/latex/images/fitted_9subject_dre.eps", device = cairo_ps, dpi = 500, width = 9, height = 8)
plotFittedDRESubject = function(modelObject, weighted = F, pid=NA,
                         FONT_SIZE=10, showTitle=T){
  data.id = modelObject$model_info$coxph_components$data
  
  if(is.na(pid)){
    pid = sample(data.id$P_ID, size = 1)
    print(paste("Choosing patient", pid, "randomly because pid not provided"))
  }
  
  data = modelObject$model_info$mvglmer_components$data
  data = data[data$P_ID == pid,]
  lastBiopsyTime =  max(data$visitTimeYears[!is.na(data$gleason)])
  data = data[!is.na(data$high_dre),] 
  
  rowNums = modelObject$model_info$mvglmer_components$id1
  dreFit = if(weighted==F){
    fitted(modelObject, process = "Longitudinal", type="Subject")[[1]]
  }else{
    fitted_weighted(modelObject)[[1]]
  }
  dreFit = dreFit[rowNums==which(data.id$P_ID==pid)]
  dreFit = plogis(dreFit)
  
  dreObserved = data$high_dre[data$P_ID==pid]
  
  plotdf = data.frame(time = data$visitTimeYears, dreObserved=dreObserved, dreFit = dreFit)
  
  plot = ggplot(data=plotdf) + 
    geom_point(aes(x=time,y=dreObserved, shape="Observed DRE"), size=3) + 
    geom_line(aes(x=time, y=dreFit, linetype="Fitted Pr (DRE > T1c)")) +
    geom_vline(aes(xintercept=lastBiopsyTime, linetype="Latest biopsy"), show.legend =  F) + 
    xlab("Follow-up time (years)") + ylab("Pr (DRE > T1c)") +
    scale_shape_manual(name="", values=18, labels="Observed DRE") +
    scale_linetype_manual(name="", values=c("dashed","solid"), 
                          labels=c("Fitted Pr (DRE > T1c)", "Latest biopsy")) +
    theme_bw() + 
    scale_y_continuous(breaks = seq(0,1, by = 0.25), labels=paste0(seq(0,1,by=0.25)*100, "%"),limits = c(0,1))+
    theme(text = element_text(size=FONT_SIZE), axis.text=element_text(size=FONT_SIZE),
          legend.text = element_text(size=FONT_SIZE),
          axis.line = element_line(),
          panel.grid.minor = element_blank()) 
  
  if(showTitle){
    plot + ggtitle(paste("Patient", pid)) 
  }else{
    plot
  }
}

#set.seed(2018)
#temp = lapply(1:9, plotFittedPSASubject, modelObject=mvJoint_dre_psa_dre_value, FONT_SIZE=13, showTitle=F)
#subjectplot = ggpubr::ggarrange(plotlist = temp, nrow = 3, ncol=3, common.legend = T, legend = "bottom")
#ggsave(subjectplot, file="report/decision_analytic/mdm/latex/images/fitted_9subject_psa.eps", device = cairo_ps, dpi = 500, width = 9, height = 8)
plotFittedPSASubject = function(modelObject, weighted=F, pid=NA,
                                FONT_SIZE=10, showTitle=T){
  data.id = modelObject$model_info$coxph_components$data
  
  if(is.na(pid)){
    pid = sample(data.id$P_ID, size = 1)
    print(paste("Choosing patient", pid, "randomly because pid not provided"))
  }
  
  data = modelObject$model_info$mvglmer_components$data
  data = data[data$P_ID == pid,]
  lastBiopsyTime =  max(data$visitTimeYears[!is.na(data$gleason)])
  data = data[!is.na(data$log2psaplus1),] 
  
  rowNums = modelObject$model_info$mvglmer_components$id2
  psaFit = if(weighted==F){
    fitted(modelObject, process = "Longitudinal", type="Subject")[[2]]
  }else{
    fitted_weighted(modelObject)[[2]]
  }
  psaFit = psaFit[rowNums==which(data.id$P_ID==pid)]
  
  log2psaplus1Observed = data$log2psaplus1[data$P_ID==pid]
  
  plotdf = data.frame(time = data$visitTimeYears, log2psaplus1Observed=log2psaplus1Observed, psaFit = psaFit)
  plot = ggplot(data=plotdf) + 
    geom_point(aes(x=time,y=log2psaplus1Observed, shape="Observed PSA"), size=2) + 
    geom_line(aes(x=time, y=psaFit, linetype="Fitted PSA")) + 
    geom_vline(aes(xintercept=lastBiopsyTime, linetype="Latest biopsy"), show.legend =  F) + 
    scale_shape_manual(name="", values=16, labels=expression('Observed log'[2]*'(PSA + 1)')) +
    scale_linetype_manual(name="", values=c("dashed","solid"), 
                          labels=c(expression('Fitted log'[2]*'(PSA + 1)'), "Latest biopsy")) +
    theme_bw() + 
    theme(text = element_text(size=FONT_SIZE), axis.text=element_text(size=FONT_SIZE),
          legend.text = element_text(size=FONT_SIZE),
          axis.line = element_line(),
          panel.grid.minor = element_blank()) +  
    xlab("Follow-up time (years)") + ylab(expression('log'[2]*'(PSA + 1)')) 
  
  if(showTitle){
    plot + ggtitle(paste("Patient", pid)) 
  }else{
    plot
  }
}

#ggsave(plotFittedPSAMarginal(mvJoint_dre_psa_dre_value, FONT_SIZE = 9), file="report/decision_analytic/mdm/latex/images/marginal_psa.eps", device = cairo_ps, dpi = 500, height=4.5/1.3333, width=4.5)
plotFittedPSAMarginal = function(modelObject, weighted=F, FONT_SIZE=10){
  visitTimeYears = seq(0,10, 0.5)
  df = data.frame(Age = 70, visitTimeYears)
  betas_psa_matrix = modelObject$mcmc$betas2
  
  fixedPSAFormula = ~ 1 +I(Age - 70) +  I((Age - 70)^2) + ns(visitTimeYears, knots=c(0.1, 0.7, 4), Boundary.knots=c(0, 5.42))
  
  xbeta = model.matrix(fixedPSAFormula, df) %*% t(betas_psa_matrix)
  
  meanXbeta = apply(xbeta, 1, mean)
  upperXbeta = apply(xbeta, 1, quantile, probs=0.975)
  lowerXbeta = apply(xbeta, 1, quantile, probs=0.025)
  
  ggplot() + 
    geom_line(aes(x=visitTimeYears, y=meanXbeta)) +
    geom_ribbon(aes(x=visitTimeYears, ymin = lowerXbeta, ymax = upperXbeta), fill="grey", alpha=0.5) + 
    theme_bw() + 
    theme(text = element_text(size=FONT_SIZE), axis.text=element_text(size=FONT_SIZE),
          axis.line = element_line(),
          panel.grid.minor = element_blank()) +  
    xlab("Follow-up time (years)") + ylab(expression('log'[2]*'(PSA + 1)'))
}

plotPSAResidualVsFitted = function(modelObject, weighted=T){
  data = modelObject$model_info$mvglmer_components$data
  data = data[!is.na(data$log2psaplus1),] 
  
  psaFit = if(weighted==F){
    fitted(modelObject, process = "Longitudinal", type="Subject")[[2]]
  }else{
    fitted_weighted(modelObject)[[2]]
  }
  log2psaplus1Observed = data$log2psaplus1
  residualPSA = log2psaplus1Observed - psaFit
  
  plot = qplot(y=residualPSA, x=log2psaplus1Observed, xlab = "Fitted", ylab = "Residuals")
  print(plot)
  return(plot)
}

#change df=100 when trying with a model which uses normality assumption on error distribution
#ggsave(qqplotPSA(mvJoint_dre_psa_dre_value, FONT_SIZE = 9), file="report/decision_analytic/mdm/latex/images/qqplot.eps", device = cairo_ps, dpi = 500, height=4.5/1.3333, width=4.5)
qqplotPSA = function(modelObject, df = 3, weighted=F, probs=c(0.25, 0.75),
                     FONT_SIZE=15){
  
  data = modelObject$model_info$mvglmer_components$data
  data = data[!is.na(data$log2psaplus1),] 
  
  psaFit = if(weighted==F){
    fitted(modelObject, process = "Longitudinal", type="Subject")[[2]]
  }else{
    fitted_weighted(modelObject)[[2]]
  }
  log2psaplus1Observed = data$log2psaplus1
  residualPSA = log2psaplus1Observed - psaFit
  
  residualPSA_quantiles <- quantile(residualPSA, probs, names = FALSE, type = 7, na.rm = TRUE)
  theoretical_quantiles = qt(probs, df=df)
  slope <- diff(residualPSA_quantiles)/diff(theoretical_quantiles)
  intercept = residualPSA_quantiles[1L] - slope * theoretical_quantiles[1L]
  
  ggplot() + geom_qq(aes(sample=residualPSA), 
                     dparams = list(df=df),
                     distribution = qt) + 
    geom_abline(intercept = intercept, slope = slope) + 
    theme_bw() + 
    theme(text = element_text(size=FONT_SIZE), axis.text=element_text(size=FONT_SIZE),
          axis.line = element_line(),
          panel.grid.minor = element_blank()) +  
    xlab("T-distribution (df=3) quantiles") + ylab("Residual quantiles")
}

getFittedPSAVelocities = function(fittedJointModel){
  generateTruePSASlope = function(visitTimeYears, randomEff_psa_slope){
    betas_psa_time = fittedJointModel$statistics$postMeans$betas2[4:7]
    
    fixedPSASlopeFormula = ~ 0 + dns(visitTimeYears, knots=c(0.1, 0.7, 4), Boundary.knots=c(0, 5.42))
    randomPSASlopeFormula = ~ 0 + dns(visitTimeYears, knots=c(0.1, 0.7, 4), Boundary.knots=c(0, 5.42))
    
    df = data.frame(visitTimeYears)
    model.matrix(fixedPSASlopeFormula, df) %*% betas_psa_time + model.matrix(randomPSASlopeFormula, df) %*% as.numeric(randomEff_psa_slope)
  }
  
  b_subjects = fittedJointModel$statistics$postMeans$b
  
  data=fittedJointModel$model_info$mvglmer_components$data
  sapply(1:length(unique(data$P_ID)), function(patientRow){
    randomEff_psa_slope = b_subjects[patientRow,4:7]
    uniquePIDs = unique(data$P_ID)
    generateTruePSASlope(data$visitTimeYears[data$P_ID==uniquePIDs[patientRow]], randomEff_psa_slope)
  })
}
