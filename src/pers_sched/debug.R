plotVarianceOverTime = function(modelObject, prias_P_ID, sd=T){
  library(latex2exp)
  source("../JMBayes/Anirudh/dev/varCondFailureTime.R")
  environment(varCondFailureTime) <- asNamespace('JMbayes')
  
  ct= makeCluster(detectCores())
  registerDoParallel(ct)
  
  prias_long_i = prias_long[prias_long$P_ID==prias_P_ID,]
  
  prias_long_i$variances=foreach(k=1:nrow(prias_long_i),.combine='c', .export=c("varCondFailureTime"),
                                 .packages=c("JMbayes", "splines")) %dopar%{
                                   #variances=foreach(k=1:5,.combine='rbind', .packages=c("JMbayes", "splines")) %dopar%{
                                   subset = prias_long_i[1:k,]
                                   last.time = tail(subset$visitTimeYears[!is.na(subset$gleason)],1)
                                   return(varCondFailureTime(modelObject, subset[!is.na(subset$psa),], "P_ID", last.time, maxPossibleFailureTime = 20))
                                 }
  
  stopCluster(ct)
  
  biopsyTimes = prias_long_i$visitTimeYears[!is.na(prias_long_i$gleason)]
  prias_long_i$biopsyTimes = c(biopsyTimes, rep(NA, nrow(prias_long_i)-length(biopsyTimes)))
  
  pp = ggplot(data=prias_long_i) + 
    geom_line(aes(x=visitTimeYears, y=sqrt(variances), linetype="Fitted")) + 
    geom_vline(aes(xintercept = 0, linetype="Exp. GR Time"), show.legend = F) +
    geom_vline(aes(xintercept = 0, linetype="Dyn. risk GR"), show.legend = F) +
    geom_vline(aes(xintercept = biopsyTimes, na.rm=T, linetype="Latest biopsy"), show.legend = F) + 
    xlab("Time (Years)") + 
    scale_linetype_manual(values=c("dotted", "twodash", "dashed", "solid"))  +
    ylab(TeX('$SD_g\\left[T^*_j\\right]$')) + ticksX(0, max = 20, 1) +
    ticksY(0, 10, 1) + ylim(0,7) +
    theme(text = element_text(size=11), axis.text=element_text(size=11), 
          legend.title = element_blank(), legend.text = element_text(size=11),
          plot.title = element_text(hjust = 0.5)) 
  
  print(pp)
  return(pp) 
}