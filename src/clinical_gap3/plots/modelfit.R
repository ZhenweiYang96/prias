library(JMbayes)
library(splines)
library(ggplot2)
library(ggpubr)

load("Rdata/gap3/PRIAS_2019/mvJoint_psa_time_scaled.Rdata")
load("Rdata/gap3/PRIAS_2019/cleandata.Rdata")

SUCCESS_COLOR = 'green'
DANGER_COLOR = 'red'
THEME_COLOR = 'dodgerblue4'
MAX_FOLLOW_UP = 10
POINT_SIZE = 2
FONT_SIZE = 13

#change df=100 when trying with a model which uses normality assumption on error distribution
qqplotPSA = function(modelObject, df = 3, probs=c(0.25, 0.75),
                     normalDist = F){
  
  if(class(modelObject)=="mvglmer"){
    data = modelObject$data
    data = data[!is.na(data$log2psaplus1),] 
    
    fittedmarginal = modelObject$components$X2 %*% modelObject$postMeans$betas2
    randEff = modelObject$postMeans$b[as.numeric(droplevels(data$P_ID)),3:7,drop=FALSE]
    fittedSubjects = rowSums(modelObject$components$Z2 * randEff)
    psaFit = fittedmarginal + fittedSubjects
    
  }else{
    data = modelObject$model_info$mvglmer_components$data
    data = data[!is.na(data$log2psaplus1),] 
    
    psaFit = fitted(modelObject, process = "Longitudinal", type="Subject")[[1]]
  }
  
  log2psaplus1Observed = data$log2psaplus1
  residualPSA = log2psaplus1Observed - psaFit
  
  residualPSA_quantiles <- quantile(residualPSA, probs, names = FALSE, type = 7, na.rm = TRUE)
  if(normalDist==T){
    theoretical_quantiles = qnorm(probs)
  }else{
    theoretical_quantiles = qt(probs, df=df)  
  }
  
  slope <- diff(residualPSA_quantiles)/diff(theoretical_quantiles)
  intercept = residualPSA_quantiles[1L] - slope * theoretical_quantiles[1L]
  
  if(normalDist==T){
    plot = ggplot() + geom_qq(aes(sample=residualPSA), 
                              distribution = qnorm) + 
      geom_abline(intercept = intercept, slope = slope) + 
      theme_bw() + 
      theme(text = element_text(size=FONT_SIZE), axis.text=element_text(size=FONT_SIZE),
            axis.line = element_line(),
            panel.grid.minor = element_blank()) +  
      xlab("Normal distribution quantiles") + ylab("Residual quantiles")
  }else{
    plot = ggplot() + geom_qq(aes(sample=residualPSA), 
                              dparams = list(df=df),
                              distribution = qt) + 
      geom_abline(intercept = intercept, slope = slope) + 
      theme_bw() + 
      theme(text = element_text(size=FONT_SIZE), axis.text=element_text(size=FONT_SIZE),
            axis.line = element_line(),
            panel.grid.minor = element_blank()) +  
      xlab("t-distribution (df=3) quantiles") + ylab("Residual quantiles")
  }
  return(plot)
}

ggsave(ggpubr::ggarrange(qqplotPSA(mvJoint_psa_time_scaled) +
                           ggtitle("Error distribution: t (df=3)"),
                         qqplotPSA(mvJoint_psa_time_scaled, normalDist = T) +
                           ggtitle("Error distribution: normal"), labels = "AUTO"),
       file="report/clinical/images/qqplot.eps",
       width = 8, height=5.5,
       device = cairo_ps, dpi = 500)

###############################################
plotFittedPSASubject = function(modelObject, pid=NA, showTitle=T){
  data.id = modelObject$model_info$coxph_components$data
  if(is.na(pid)){
    repeat{
      pid = sample(data.id$P_ID, size = 1)
      print(paste("Choosing patient", pid, "randomly because pid not provided"))
      
      data = modelObject$model_info$mvglmer_components$data
      data = data[data$P_ID == pid,]
      lastBiopsyTime = max(data$year_visit[!is.na(data$gleason_sum)])
      data = data[!is.na(data$log2psaplus1),]
      if(nrow(data)>4){
        break
      }
    }
  }
  
  rowNums = modelObject$model_info$mvglmer_components$id1
  psaFit = fitted(modelObject, process = "Longitudinal", type="Subject")[[1]]
  psaFit = psaFit[rowNums==which(data.id$P_ID==pid)]
  
  log2psaplus1Observed = data$log2psaplus1[data$P_ID==pid]
  
  plotdf = data.frame(time = data$year_visit, log2psaplus1Observed=log2psaplus1Observed, psaFit = psaFit)
  
  yrange = range(c(plotdf[,-1]))
  
  plot = ggplot(data=plotdf) + 
    geom_point(aes(x=time,y=log2psaplus1Observed, shape="Observed PSA"), size=POINT_SIZE, color=THEME_COLOR) + 
    geom_line(aes(x=time, y=psaFit, linetype="Fitted PSA"), color=THEME_COLOR) + 
    geom_vline(aes(xintercept=lastBiopsyTime, linetype="Latest biopsy"), show.legend =  F) + 
    scale_shape_manual(name="", values=16, labels=expression('Observed log'[2]*'(PSA + 1)')) +
    scale_linetype_manual(name="", values=c("dashed","solid"), 
                          labels=c(expression('Fitted log'[2]*'(PSA + 1)'), "Latest biopsy")) +
    theme_bw() + 
    theme(text = element_text(size=FONT_SIZE), 
          axis.text=element_text(size=FONT_SIZE),
          axis.text.y = element_text(size=FONT_SIZE, color=THEME_COLOR),
          axis.title.y = element_text(size=FONT_SIZE, color=THEME_COLOR),
          axis.line = element_line(),
          panel.grid.minor = element_blank()) +  
    scale_y_continuous(breaks = seq(yrange[1], yrange[2], length.out = 4),
                       labels = round(seq(yrange[1], yrange[2], length.out = 4),1),
                       limits = yrange)+
    xlab("Follow-up time (years)") + ylab(expression('log'[2]*'(PSA + 1)')) 
  
  if(showTitle){
    plot + ggtitle(paste("Patient", pid)) 
  }else{
    plot
  }
}

set.seed(2019)
temp = list()
for(m in 1:9){
  temp[[m]] = plotFittedPSASubject(modelObject=mvJoint_psa_time_scaled, pid=NA, showTitle=F)
}
subjectplot = ggpubr::ggarrange(plotlist = temp, nrow = 3, ncol=3, common.legend = T, legend = "bottom")
ggsave(subjectplot, file="report/clinical/images/fitted_9subject_psa.eps", device = cairo_ps, 
       dpi = 500, width = 8, height = 8)
