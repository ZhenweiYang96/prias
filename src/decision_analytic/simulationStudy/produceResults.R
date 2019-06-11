#First load all the results and extract the actual schedule results
resFiles = list.files("/home/a_tomer/Data/final_res_2nd_paper/", full.names = T)

runFixedSchedule = function(progressionTimes, biopsyTimes){
  if(!10 %in% biopsyTimes){
    biopsyTimes = c(biopsyTimes, 10)
  }
  
  nb = sapply(progressionTimes, function(x){which(biopsyTimes>=x)[1]})
  detectionTimes = biopsyTimes[nb]
  offset = detectionTimes - progressionTimes
  
  return(cbind(nb, offset))
}

library(doParallel)
simResults_n500_method_all = foreach(i=1:50, .combine="rbind") %do%{
  load(resFiles[i])
  
  print(i)
  #print(resFiles[i])
  
  #seed = jointModelData$seed
  #set.seed(seed)
  
  #real_f1 = jointModelData$testData$scheduleResults[jointModelData$testData$scheduleResults$methodName==KAPPAF1Score,]
  #real_f1$methodName = "Risk (F1) Real"
  
  #real_f1[,c("nb", "offset")] = runDynRiskGRSchedule(jointModelData, "F1score")
  #jointModelData$testData$scheduleResults = rbind(jointModelData$testData$scheduleResults, real_f1)
  
  #save(jointModelData, file=resFiles[i])
  
  #month_12=runFixedSchedule(jointModelData$testData$testDs.id$progression_time, seq(1, 10, 1))
  month_18=runFixedSchedule(jointModelData$testData$testDs.id$progression_time, seq(1.5, 10, 1.5))
  month_24=runFixedSchedule(jointModelData$testData$testDs.id$progression_time, seq(2, 10, 2))
  month_24s=runFixedSchedule(jointModelData$testData$testDs.id$progression_time, c(1, 10, 2))
  #month_30=runFixedSchedule(jointModelData$testData$testDs.id$progression_time, seq(2.5, 10, 2.5))
  month_36=runFixedSchedule(jointModelData$testData$testDs.id$progression_time, seq(3, 10, 3))
  #month_42=runFixedSchedule(jointModelData$testData$testDs.id$progression_time, seq(3.5, 10, 3.5))
  month_48=runFixedSchedule(jointModelData$testData$testDs.id$progression_time, seq(4, 10, 4))
  
  res = jointModelData$testData$scheduleResults
  res$nb[res$methodName=="Biennial"] = month_24[,1]
  res$offset[res$methodName=="Biennial"] = month_24[,2]
  
  res$nb[res$methodName=="18 Months"] = month_18[,1]
  res$offset[res$methodName=="18 Months"] = month_18[,2]
  
  res_24s = res[res$methodName=="Biennial",]
  res_24s$methodName = "Biennial_S"
  res_24s$nb = month_24s[,1]
  res_24s$offset = month_24s[,2]
  
  res_36 = res[res$methodName=="Biennial",]
  res_36$methodName = "Triennial"
  res_36$nb = month_36[,1]
  res_36$offset = month_36[,2]
  
  res_48 = res[res$methodName=="Biennial",]
  res_48$methodName = "Quadriennial"
  res_48$nb = month_48[,1]
  res_48$offset = month_48[,2]
  
  res = rbind(res, res_24s, res_36, res_48)

  ret = res
  rm(jointModelData)
  return(ret)
}

scheduleResCombined = simResults_n500_method_all
scheduleResCombined = simResults_n500_method_all[simResults_n500_method_all$methodName %in% c("Annual", "Risk (10%)", "Risk (5%)", "Risk (F1) Real", "PRIAS"),]
scheduleResCombined = simResults_n500_method_all[simResults_n500_method_all$methodName %in% c("Annual", "18 Months", "Biennial", "Triennial", "Quadriennial",
                                                                                              "PRIAS", "Risk (5%)", "Risk (10%)", "Risk (15%)", "Risk (F1) Real"),]
scheduleResCombined$methodName = droplevels(scheduleResCombined$methodName)

scheduleResCombined$nb[scheduleResCombined$offset < 0] = scheduleResCombined$nb[scheduleResCombined$offset<0] + 1
scheduleResCombined$offset[scheduleResCombined$offset < 0] = 10 - scheduleResCombined$progression_time[scheduleResCombined$offset<0]

levels(scheduleResCombined$methodName)[3] = "Risk: 10%"
levels(scheduleResCombined$methodName)[4] = "Risk: 5%"
levels(scheduleResCombined$methodName)[5] = "Risk: F1"


levels(scheduleResCombined$methodName)[5] = "Risk: 10%"
levels(scheduleResCombined$methodName)[6] = "Risk: 15%"
levels(scheduleResCombined$methodName)[7] = "Risk: 5%"
levels(scheduleResCombined$methodName)[8] = "Risk: F1"

FONT_SIZE = 11

getBoxplotStatsDf=function(progression_time_low, progression_time_high, attribute){
  temp = scheduleResCombined[scheduleResCombined$progression_time>=progression_time_low & scheduleResCombined$progression_time<=progression_time_high,]
  res = unlist(by(temp$methodName, data = temp[, attribute], FUN = function(x){
    boxplot.stats(x)$stats
  }))
  
  resDf = data.frame(matrix(res, ncol=5, byrow = T))
  resDf$methodName = levels(scheduleResCombined$methodName)
  return(resDf)
}

temp = getBoxplotStatsDf(0,10, "nb")[,c(6,2,3,4)]
temp = cbind(temp, round(getBoxplotStatsDf(0, 9.999999, "offset")[, 2:4], 1))
write.csv(temp, file="all.csv", row.names = F)

temp = getBoxplotStatsDf(0,3.5, "nb")[,c(6,2,3,4)]
temp = cbind(temp, round(getBoxplotStatsDf(0, 3.5, "offset")[, 2:4], 1))
write.csv(temp, file="fast.csv", row.names = F)

temp = getBoxplotStatsDf(3.50001,9.999999, "nb")[,c(6,2,3,4)]
temp = cbind(temp, round(getBoxplotStatsDf(3.50001,9.999999, "offset")[, 2:4], 1))
write.csv(temp, file="intermediate.csv", row.names = F)

temp = getBoxplotStatsDf(10,10, "nb")[,c(6,2,3,4)]
temp = cbind(temp, round(getBoxplotStatsDf(10,10, "offset")[, 2:4], 1))
write.csv(temp, file="slow.csv", row.names = F)

MEDIAN_WIDTH_BOXPLOT = 2

gfastNb = ggplot(data=getBoxplotStatsDf(0,3.5, "nb")) + 
  geom_boxplot(aes(ymin = X1, lower = X2, middle = X3, upper = X4, ymax = X5, x=methodName),
               stat = "identity", fatten=3) + coord_flip() + 
  scale_y_continuous(breaks=seq(1,10, length.out = 4), limits = c(1,10)) +
  theme_bw() + 
  theme(text = element_text(size=FONT_SIZE), axis.text=element_text(size=FONT_SIZE),
        axis.line = element_line(),
        axis.text.y = element_text(size=FONT_SIZE, color = "gray30"),
        axis.title.y = element_text(size=FONT_SIZE, color = "black"),
        axis.text.x = element_text(size=FONT_SIZE, color="gray40"),
        legend.background = element_blank(), legend.position = "top",
        legend.text = element_text(size=FONT_SIZE-3), title = element_text(size=FONT_SIZE-1))  +
  xlab("Schedule") + ylab("Number of Biopsies") + 
  ggtitle("Fast progressing (30% patients)")

gfastOffset = ggplot(data=getBoxplotStatsDf(0,3.5, "offset")) + 
  geom_boxplot(aes(ymin = X1, lower = X2, middle = X3, upper = X4, ymax = X5, x=methodName),
               stat = "identity", fatten=MEDIAN_WIDTH_BOXPLOT) + coord_flip() + 
  scale_y_continuous(breaks=0:5, limits = c(0,5.5)) +
  theme_bw() + 
  theme(text = element_text(size=FONT_SIZE), axis.text=element_text(size=FONT_SIZE),
        axis.line = element_line(),
        #axis.text.y = element_text(size=FONT_SIZE, color = "gray30"),
        axis.text.y = element_blank(),
        #axis.title.y = element_text(size=FONT_SIZE, color = "black"),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size=FONT_SIZE, color="gray40"),
        legend.background = element_blank(), legend.position = "top",
        legend.text = element_text(size=FONT_SIZE-3)) + 
  xlab("Schedule") + ylab("Delay in detection of\ncancer progression (years)")
#  geom_hline(yintercept = 2, linetype="dashed")

gFast = ggpubr::ggarrange(gfastNb, gfastOffset, ncol=2, widths = c(1.4,1), align="h")

gIntermediateNb = ggplot(data=getBoxplotStatsDf(3.50001,9.999999, "nb")) + 
  geom_boxplot(aes(ymin = X1, lower = X2, middle = X3, upper = X4, ymax = X5, x=methodName),
               stat = "identity", fatten=MEDIAN_WIDTH_BOXPLOT) + coord_flip() + 
  scale_y_continuous(breaks=seq(1,10, length.out = 4)) +
  theme_bw() + 
  theme(text = element_text(size=FONT_SIZE), axis.text=element_text(size=FONT_SIZE),
        axis.line = element_line(),
        axis.text.y = element_text(size=FONT_SIZE, color = "gray30"),
        axis.title.y = element_text(size=FONT_SIZE, color = "black"),
        axis.text.x = element_text(size=FONT_SIZE, color="gray40"),
        legend.background = element_blank(), legend.position = "top",
        legend.text = element_text(size=FONT_SIZE-3), title = element_text(size=FONT_SIZE-1))  +
  xlab("Schedule") + ylab("Number of Biopsies") + 
  ggtitle("Intermediate progressing (20% patients)")

gIntermediateOffset = ggplot(data=getBoxplotStatsDf(3.50001,9.999999, "offset")) + 
  geom_boxplot(aes(ymin = X1, lower = X2, middle = X3, upper = X4, ymax = X5, x=methodName),
               stat = "identity", fatten=MEDIAN_WIDTH_BOXPLOT) + coord_flip() + 
  scale_y_continuous(breaks=0:5, limits = c(0,5.5)) +
  theme_bw() + 
  theme(text = element_text(size=FONT_SIZE), axis.text=element_text(size=FONT_SIZE),
        axis.line = element_line(),
        #axis.text.y = element_text(size=FONT_SIZE, color = "gray30"),
        axis.text.y = element_blank(),
        #axis.title.y = element_text(size=FONT_SIZE, color = "black"),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size=FONT_SIZE, color="gray40"),
        legend.background = element_blank(), legend.position = "top",
        legend.text = element_text(size=FONT_SIZE-3)) + 
  xlab("Schedule") + ylab("Delay in detection of\ncancer progression (years)")
#  geom_hline(yintercept = 2, linetype="dashed")

gIntermediate = ggpubr::ggarrange(gIntermediateNb, gIntermediateOffset, ncol=2, widths = c(1.4,1), align = "h")

gSlowNb = ggplot(data=getBoxplotStatsDf(10,10, "nb")) + 
  geom_boxplot(aes(ymin = X1, lower = X2, middle = X3, upper = X4, ymax = X5, x=methodName),
               stat = "identity", fatten=MEDIAN_WIDTH_BOXPLOT) + coord_flip() + 
  scale_y_continuous(breaks=c(1,4,7,10)) +
  theme_bw() + 
  theme(text = element_text(size=FONT_SIZE), axis.text=element_text(size=FONT_SIZE),
        axis.line = element_line(),
        axis.text.y = element_text(size=FONT_SIZE, color = "gray30"),
        axis.title.y = element_text(size=FONT_SIZE, color = "black"),
        axis.text.x = element_text(size=FONT_SIZE, color="gray40"),
        legend.background = element_blank(), legend.position = "top",
        legend.text = element_text(size=FONT_SIZE-3),title = element_text(size=FONT_SIZE-1))  +
  xlab("Schedule") + ylab("Number of Biopsies") + 
  ggtitle("Slow progressing (50% patients)")

gSlowOffset = ggplot() + geom_text(aes(x=1, y=2, 
                                       label="No cancer progression\nobserved in 10 year follow-up period\nfor patients with slow speed of cancer\nprogression. Hence, no boxplot for\ndelay in detection of cancer\nprogression.")) + theme_bw() + 
  theme(axis.ticks = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank())

gSlow = ggpubr::ggarrange(gSlowNb, gSlowOffset, ncol=2, widths = c(1.4,1), align="h")

combinedPlot = ggpubr::ggarrange(gFast, gIntermediate,gSlow, nrow=3, ncol=1, 
                  labels = "AUTO")
ggsave(combinedPlot, filename = "report/decision_analytic/mdm/latex/images/sim_res_combined.eps",
       height=8, width=7, dpi = 500)


#How F1 score works
patientDs = prias_long[prias_long$P_ID==2340 & prias_long$visitTimeYears<=4.1,]
f1score=ggplot() + geom_point(data=patientDs,aes(x=visitTimeYears, y=log2psa),size=2, 
                      color='dodgerblue4') +
  geom_vline(xintercept = 2.6) + 
  geom_segment(aes(x=4, xend=4, y=-Inf, yend=4.15)) +
  geom_ribbon(data=data.frame(visitTimeYears=seq(2.6, 4, 0.1)), 
              aes(ymin=2.722466, ymax=4.5, x=visitTimeYears, fill="firebrick2", alpha=0.5)) +
  xlab("Follow-up time (years)") + ylab(expression('log'[2]*'(PSA)')) +
  theme(text = element_text(size=15), axis.text=element_text(size=15),
        legend.position = "none",
        axis.text.y = element_text(size=15, color = "dodgerblue4"),
        axis.title.y = element_text(size=15, color = "dodgerblue4"),
        axis.title.y.right = element_text(size=15, color = "firebrick2"),
        axis.text.y.right = element_text(size=15, color = "firebrick2")) +
  scale_x_continuous(breaks = c(0, 1.5, 2.6, 4, 6.5), limits = c(0, 6.5), 
                     labels = str_wrap(c("0", "1.5", "t=2.6 (Latest biopsy)", "s=4 (Current Visit)", "6.5"), width=8))
  
ggsave(.Last.value, file="f1score(s-t).eps", width = 8.27, device = cairo_ps)
