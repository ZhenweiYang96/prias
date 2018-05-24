source("src/decision_analytic/fittingModel/weighted_fitted.R")

plotFittedDRE = function(modelObject, weighted = T, pid=NA){
  if(is.na(pid)){
    pid = sample(prias.id$P_ID, size = 1)
    print(paste("Choosing patient", pid, "randomly because pid not provided"))
  }
  
  data = modelObject$model_info$mvglmer_components$data
  data = data[!is.na(data$high_dre),] 
  
  rowNums = modelObject$model_info$mvglmer_components$id1
  dreFit = if(weighted==F){
    fitted(modelObject, process = "Longitudinal", type="Subject")[[1]]
  }else{
    fitted_weighted(modelObject)[[1]]
  }
  dreFit = dreFit[rowNums==which(prias.id$P_ID==pid)]
  dreFit = 1/(1 + exp(-(dreFit)))
  
  dreObserved = data$high_dre[data$P_ID==pid]
  
  plotdf = data.frame(time = data$visitTimeYears[data$P_ID == pid], dreObserved=dreObserved, dreFit = dreFit)
  ggplot(data=plotdf) + geom_point(aes(x=time,y=dreObserved)) + geom_line(aes(x=time, y=dreFit)) +
    ylim(c(0,1)) + xlab("Time (years)") + ylab("DRE > T1c") + ggtitle(paste("Patient", pid)) 
  
}

plotFittedPSA = function(modelObject, weighted=T, pid=NA){
  if(is.na(pid)){
    pid = sample(prias.id$P_ID, size = 1)
    print(paste("Choosing patient", pid, "randomly because pid not provided"))
  }
  
  data = modelObject$model_info$mvglmer_components$data
  data = data[!is.na(data$log2psaplus1),] 
  
  rowNums = modelObject$model_info$mvglmer_components$id2
  psaFit = if(weighted==F){
    fitted(modelObject, process = "Longitudinal", type="Subject")[[2]]
  }else{
    fitted_weighted(modelObject)[[2]]
  }
  psaFit = psaFit[rowNums==which(prias.id$P_ID==pid)]
  
  log2psaplus1Observed = data$log2psaplus1[data$P_ID==pid]
  
  plotdf = data.frame(time = data$visitTimeYears[data$P_ID == pid], log2psaplus1Observed=log2psaplus1Observed, psaFit = psaFit)
  ggplot(data=plotdf) + geom_point(aes(x=time,y=log2psaplus1Observed)) + geom_line(aes(x=time, y=psaFit)) +
    xlab("Time (years)") + ylab(expression('log'[2]*' (PSA + 1)')) + ggtitle(paste("Patient", pid)) 
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
qqplotPSA = function(modelObject, df = 3, weighted=T, probs=c(0.25, 0.75)){
  
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
    xlab("T-distribution (df=3) quantiles") + ylab("Residual quantiles") + 
    theme(text = element_text(size=11), axis.text=element_text(size=11), 
          plot.title = element_text(hjust = 0.5, size=13)) 
}

