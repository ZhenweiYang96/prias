compareVarianceOverTime = function(jmfitNorm, jmfitT3, prias_P_ID, sd=T){
  source("../JMBayes/Anirudh/dev/varCondFailureTime.R")
  environment(varCondFailureTime) <- asNamespace('JMbayes')
  
  ct= makeCluster(detectCores())
  registerDoParallel(ct)
  
  prias_long_i = prias_long[prias_long$P_ID==prias_P_ID,]
  biopsyTimes = prias_long_i$visitTimeYears[!is.na(prias_long_i$gleason)]
  prias_long_i$biopsyTimes = c(biopsyTimes, rep(NA, nrow(prias_long_i)-length(biopsyTimes)))
  
  nrowOrig = nrow(prias_long_i)
  
  variances_norm=foreach(k=1:nrow(prias_long_i),.combine='c', .export=c("varCondFailureTime"),
                         .packages=c("JMbayes", "splines")) %dopar%{
                           #variances=foreach(k=1:5,.combine='rbind', .packages=c("JMbayes", "splines")) %dopar%{
                           subset = prias_long_i[1:k,]
                           last.time = tail(subset$visitTimeYears[!is.na(subset$gleason)],1)
                           return(varCondFailureTime(jmfitNorm, subset[!is.na(subset$psa),], "P_ID", last.time, maxPossibleFailureTime = 20))
                         }
  
  variances_t3=foreach(k=1:nrow(prias_long_i),.combine='c', .export=c("varCondFailureTime"),
                       .packages=c("JMbayes", "splines")) %dopar%{
                         #variances=foreach(k=1:5,.combine='rbind', .packages=c("JMbayes", "splines")) %dopar%{
                         subset = prias_long_i[1:k,]
                         last.time = tail(subset$visitTimeYears[!is.na(subset$gleason)],1)
                         return(varCondFailureTime(jmfitT3, subset[!is.na(subset$psa),], "P_ID", last.time, maxPossibleFailureTime = 20))
                       }
  
  stopCluster(ct)
  
  plotDf = data.frame(sd = sqrt(c(variances_norm, variances_t3)), 
                      type=rep(c("Normally distributed errors", "t-distributed (df=3) errors"),each=nrowOrig), 
                      visitTimeYears=rep(prias_long_i$visitTimeYears,2))
  plotDf$biopsyTimes = c(biopsyTimes, rep(NA, nrow(plotDf)-length(biopsyTimes)))
  
  pp = ggplot(data=plotDf) +
    geom_line(aes(x=visitTimeYears, y=sd, linetype=type)) + 
    geom_vline(aes(xintercept = biopsyTimes, na.rm=T, linetype="Latest biopsy"), show.legend = F) + 
    xlab("Time (years)") + ylab(TeX('$SD_g\\left[T^*_j\\right]$')) + ticksX(0, max = 20, 1) +
    ticksY(0, 10, 1) + 
    scale_linetype_manual(values=c("solid", "dotted", "dashed"))  +
    theme(text = element_text(size=11), axis.text=element_text(size=11), 
          legend.title = element_blank(), legend.text = element_text(size=11),
          plot.title = element_text(hjust = 0.5, size=11)) 
  
  print(pp)
  return(pp) 
}