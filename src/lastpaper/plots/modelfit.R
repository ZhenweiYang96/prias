library(JMbayes)
library(splines)
library(ggplot2)
library(ggpubr)

load("Rdata/lastpaper/fitted_model/mvJoint_dre_psa_2knots_quad_age.Rdata")
load("Rdata/gap3/PRIAS_2019/cleandata.Rdata")

SUCCESS_COLOR = 'green'
DANGER_COLOR = 'red'
THEME_COLOR = 'dodgerblue4'
MAX_FOLLOW_UP = 10
POINT_SIZE = 2
FONT_SIZE = 12

#change df=100 when trying with a model which uses normality assumption on error distribution
qqplotPSA = function(modelObject, df = 3, probs=c(0.25, 0.75),
                     normalDist = F){
  
  data = modelObject$model_info$mvglmer_components$data
  data = data[!is.na(data$log2psaplus1),] 
  
  psaFit = fitted(modelObject, process = "Longitudinal", type="Subject")[[2]]
  
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
      xlab(paste0("t-distribution (df=", df, ") quantiles")) + ylab("Residual quantiles")
  }
  return(plot)
}

ggsave(ggpubr::ggarrange(qqplotPSA(mvJoint_dre_psa_2knots_quad_age) +
                           ggtitle("t (df=3) distribution"),
                         qqplotPSA(mvJoint_dre_psa_2knots_quad_age, df = 4) +
                           ggtitle("t (df=4) distribution"),
                         qqplotPSA(mvJoint_dre_psa_2knots_quad_age, normalDist = T) +
                           ggtitle("Normal distribution"), 
                         ncol = 3, nrow = 1, labels = "AUTO"),
       file="report/lastpaper/images/qqplot.eps",
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
  
  rowNums = modelObject$model_info$mvglmer_components$id2
  psaFit = fitted(modelObject, process = "Longitudinal", type="Subject")[[2]]
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
  temp[[m]] = plotFittedPSASubject(modelObject=mvJoint_dre_psa_2knots_quad_age, pid=NA, showTitle=F)
}
subjectplot = ggpubr::ggarrange(plotlist = temp, nrow = 3, ncol=3, common.legend = T, legend = "bottom")
ggsave(subjectplot, file="report/lastpaper/images/fitted_9subject_psa.eps", device = cairo_ps, 
       dpi = 500, width = 8, height = 8)


############
plotFittedDRESubject = function(modelObject, pid=NA, showTitle=T){
  data.id = modelObject$model_info$coxph_components$data
  
  data.id = modelObject$model_info$coxph_components$data
  if(is.na(pid)){
    repeat{
      pid = sample(data.id$P_ID, size = 1)
      print(paste("Choosing patient", pid, "randomly because pid not provided"))
      
      data = modelObject$model_info$mvglmer_components$data
      data = data[data$P_ID == pid,]
      lastBiopsyTime = max(data$year_visit[!is.na(data$gleason_sum)])
      data = data[!is.na(data$palpable_dre),]
      if(nrow(data)>2){
        break
      }
    }
  }
  
  rowNums = modelObject$model_info$mvglmer_components$id1
  dreFit = fitted(modelObject, process = "Longitudinal", type="Subject")[[1]]
  dreFit = dreFit[rowNums==which(data.id$P_ID==pid)]
  dreFit = plogis(dreFit)
  
  dreObserved = as.numeric(data$palpable_dre[data$P_ID==pid])
  
  plotdf = data.frame(time = data$year_visit, dreObserved=dreObserved, dreFit = dreFit)
  
  plot = ggplot(data=plotdf) + 
    geom_point(aes(x=time,y=dreObserved, shape="Observed DRE"), size=3) + 
    geom_line(aes(x=time, y=dreFit, linetype="Fitted Pr (DRE > T1c)")) +
    geom_vline(aes(xintercept=lastBiopsyTime, linetype="Latest biopsy"), show.legend =  F) + 
    xlab("Follow-up time (years)") + ylab("Pr (DRE > T1c)") +
    scale_shape_manual(name="", values=17, labels="Observed DRE") +
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

set.seed(2019)
temp = list()
for(m in 1:9){
  temp[[m]] = plotFittedDRESubject(modelObject=mvJoint_dre_psa_2knots_quad_age, pid=NA,  showTitle=F)
}
subjectplot = ggpubr::ggarrange(plotlist = temp, nrow = 3, ncol=3, common.legend = T, legend = "bottom")
ggsave(subjectplot, file="report/lastpaper/images/fitted_9subject_dre.eps",
       device = cairo_ps, dpi = 500, width = 8, height = 8)

