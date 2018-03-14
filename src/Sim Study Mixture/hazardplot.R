library(JMbayes)
library(splines)
library(ggplot2)

#number of subjects
weibullScales = c(4,5,6)
weibullShapes = c(1.5,3,4.5)

times = seq(0, 10, length.out = 500)
hazard_g1 = (weibullShapes[1]/weibullScales[1])*(times/weibullScales[1])^(weibullShapes[1]-1)
hazard_g2 = (weibullShapes[2]/weibullScales[2])*(times/weibullScales[2])^(weibullShapes[2]-1)
hazard_g3 = (weibullShapes[3]/weibullScales[3])*(times/weibullScales[3])^(weibullShapes[3]-1)

weights = t(sapply(times, function(t){
  s1 = exp((t/weibullScales[1])^weibullShapes[1])
  s2 = exp((t/weibullScales[2])^weibullShapes[2])
  s3 = exp((t/weibullScales[3])^weibullShapes[3])
  
  return(c(s1,s2,s3)/(s1+s2+s3))
}))

theoreticalHazardMatrix = cbind(hazard_g1, hazard_g2, hazard_g3)
mixtureTheoreticalHazard = rowSums(theoreticalHazardMatrix * weights)

datasetNumbers = c(1:500)
hazardMatrix = matrix(ncol=length(times), nrow=length(datasetNumbers))

for(i in datasetNumbers){
  load(paste("C:/Users/838035/Old Sim Results/Rdata/Gleason as event/Sim Study/sc_mixed_sh_mixed_f1mixed/Dt_1/simDs",i,".Rdata", sep=""))
  temp = temp[[1]]$models$mvJoint_psa_tdboth_training

  print(paste("******** Started working on Data Set: ", i, "*******"))
  
  hazardMatrix[i, ] = splineDesign(temp$control$knots, times, 
                                   ord = temp$control$ordSpline, 
                                   outer.ok = T) %*% temp$statistics$postMeans$Bs_gammas
  
  print(paste("******** End working on Data Set: ", i, "*******"))
  
  #Save RAM
  rm(temp)
}

slope_estimate = rep(NA, 500)

for(i in datasetNumbers){
  load(paste("C:/Users/838035/Old Sim Results/Rdata/Gleason as event/Sim Study/sc_mixed_sh_mixed_f1mixed/Dt_1/simDs",i,".Rdata", sep=""))
  temp = temp[[1]]$models$mvJoint_psa_tdboth_training
  
  print(paste("******** Started working on Data Set: ", i, "*******"))
  
  slope_estimate[i] = temp$statistics$postMeans$alphas[2]
  
  print(paste("******** End working on Data Set: ", i, "*******"))
  
  #Save RAM
  rm(temp)
}

fittedHazardMean = (apply((hazardMatrix), 2, mean, na.rm=T))
fittedHazardLower = (apply((hazardMatrix), 2, quantile, probs=c(0.025), na.rm=T))
fittedHazardUpper = (apply((hazardMatrix), 2, quantile, probs=c(0.975), na.rm=T))

df = data.frame(hazard = c(log(mixtureTheoreticalHazard), fittedHazardMean),
                "Baseline Hazard"=rep(c("Theoretical log (baseline hazard)", 
                                        "Fitted log(baseline hazard)"), each=500),
                confIntLow=rep(fittedHazardLower,2),
                confIntHigh=rep(fittedHazardUpper,2),
                time = rep(times,2))

ggplot(data=df[df$time > 0.1 & df$time<=8,]) + 
  geom_line(aes(x=time, y=(hazard), linetype=Baseline.Hazard)) + 
  geom_ribbon(aes(x=time, ymin=(confIntLow), ymax=(confIntHigh)), fill = "grey", alpha=0.4) + 
  scale_x_continuous(breaks = c(0.1, seq(1, 10, 1))) + 
  scale_linetype_manual(values=c("twodash", "solid")) +
  theme(legend.title=element_blank(), legend.position = "top",
        text = element_text(size=11), axis.text=element_text(size=11),
        legend.text = element_text(size=11)) + 
  xlab("Time (years)") + ylab("log (baseline hazard)")

ggsave(file="report/pers_schedule/biometrics_submission/images/sim_study/baseline_hazard.eps", 
       width=8.27*0.75, height=8.27*0.75, device=cairo_ps)
