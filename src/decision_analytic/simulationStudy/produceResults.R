#First load all the results and extract the actual schedule results
library(doParallel)

resFiles = list.files("/home/a_tomer/Results/final_res_2nd_paper/", full.names = T)

cl = makeCluster(4)
registerDoParallel(cl)

sim_study_res = foreach(i=1:length(resFiles), .combine="rbind") %dopar%{
  load(resFiles[i])
  
  #print(resFiles[i])
  
  #seed = jointModelData$seed
  #set.seed(seed)
  
  #real_f1 = jointModelData$testData$scheduleResults[jointModelData$testData$scheduleResults$methodName==KAPPAF1Score,]
  #real_f1$methodName = "Risk (F1) Real"
  
  #real_f1[,c("nb", "offset")] = runDynRiskGRSchedule(jointModelData, "F1score")
  #jointModelData$testData$scheduleResults = rbind(jointModelData$testData$scheduleResults, real_f1)
  
  #save(jointModelData, file=resFiles[i])
  ret = jointModelData$testData$scheduleResults
  rm(jointModelData)
  return(ret)
}
stopCluster(cl)

FONT_SIZE = 11

scheduleResCombined = sim_study_res[sim_study_res$methodName %in% c(ANNUAL, PRIAS, KAPPApt85, KAPPApt95, KAPPAF1Score, "Risk (F1) Real"),]
scheduleResCombined$methodName = droplevels(scheduleResCombined$methodName)

scheduleResCombined$nb[scheduleResCombined$offset < 0] = scheduleResCombined$nb[scheduleResCombined$offset<0] + 1
scheduleResCombined$offset[scheduleResCombined$offset < 0] = 10 - scheduleResCombined$progression_time[scheduleResCombined$offset<0]

levels(scheduleResCombined$methodName)[3] = "Risk: 15%"
levels(scheduleResCombined$methodName)[4] = "Risk: 5%"
levels(scheduleResCombined$methodName)[5] = "Risk: F1 (Dt=0.5)"
levels(scheduleResCombined$methodName)[6] = "Risk: F1 (Dt=s-t)"

getBoxplotStatsDf=function(progression_time_low, progression_time_high, attribute){
  temp = scheduleResCombined[scheduleResCombined$progression_time>=progression_time_low & scheduleResCombined$progression_time<=progression_time_high,]
  res = unlist(by(temp$methodName, data = temp[, attribute], FUN = function(x){
    boxplot.stats(x)$stats
  }))
  
  resDf = data.frame(matrix(res, ncol=5, byrow = T))
  resDf$methodName = levels(scheduleResCombined$methodName)
  return(resDf)
}

gfastNb = ggplot(data=getBoxplotStatsDf(0,3.5, "nb")) + 
  geom_boxplot(aes(ymin = X1, lower = X2, middle = X3, upper = X4, ymax = X5, x=methodName),
               stat = "identity") + coord_flip() + 
  scale_y_continuous(breaks=1:3) +
  theme_bw() + 
  theme(text = element_text(size=FONT_SIZE), axis.text=element_text(size=FONT_SIZE),
        axis.line = element_line(),
        axis.text.y = element_text(size=FONT_SIZE, color = "gray30"),
        axis.title.y = element_text(size=FONT_SIZE, color = "black"),
        axis.text.x = element_text(size=FONT_SIZE, color="gray40"),
        legend.background = element_blank(), legend.position = "top",
        legend.text = element_text(size=FONT_SIZE-3))  +
  xlab("Schedule") + ylab("Number of Biopsies") + ggtitle("Progression speed: Fast")

gfastOffset = ggplot(data=getBoxplotStatsDf(0,3.5, "offset")) + 
  geom_boxplot(aes(ymin = X1, lower = X2, middle = X3, upper = X4, ymax = X5, x=methodName),
               stat = "identity") + coord_flip() + 
  scale_y_continuous(breaks=0:4, limits = c(0,4.5)) +
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
  xlab("Schedule") + ylab(str_wrap("Delay in detection of cancer progression (years)", width = 30)) +
  geom_hline(yintercept = 2, linetype="dashed")

gFast = ggpubr::ggarrange(gfastNb, gfastOffset, ncol=2, widths = c(1.4,1), align="h")

gIntermediateNb = ggplot(data=getBoxplotStatsDf(3.50001,9.999999, "nb")) + 
  geom_boxplot(aes(ymin = X1, lower = X2, middle = X3, upper = X4, ymax = X5, x=methodName),
               stat = "identity") + coord_flip() + 
  scale_y_continuous(breaks=seq(1,10, length.out = 4)) +
  theme_bw() + 
  theme(text = element_text(size=FONT_SIZE), axis.text=element_text(size=FONT_SIZE),
        axis.line = element_line(),
        axis.text.y = element_text(size=FONT_SIZE, color = "gray30"),
        axis.title.y = element_text(size=FONT_SIZE, color = "black"),
        axis.text.x = element_text(size=FONT_SIZE, color="gray40"),
        legend.background = element_blank(), legend.position = "top",
        legend.text = element_text(size=FONT_SIZE-3))  +
  xlab("Schedule") + ylab("Number of Biopsies") + ggtitle("Progression speed: Intermediate")

gIntermediateOffset = ggplot(data=getBoxplotStatsDf(3.50001,9.999999, "offset")) + 
  geom_boxplot(aes(ymin = X1, lower = X2, middle = X3, upper = X4, ymax = X5, x=methodName),
               stat = "identity") + coord_flip() + 
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
  xlab("Schedule") + ylab(str_wrap("Delay in detection of cancer progression (years)", width = 30)) +
  geom_hline(yintercept = 2, linetype="dashed")

gIntermediate = ggpubr::ggarrange(gIntermediateNb, gIntermediateOffset, ncol=2, widths = c(1.4,1), align = "h")

gSlowNb = ggplot(data=getBoxplotStatsDf(10,10, "nb")) + 
  geom_boxplot(aes(ymin = X1, lower = X2, middle = X3, upper = X4, ymax = X5, x=methodName),
               stat = "identity") + coord_flip() + 
  scale_y_continuous(breaks=c(1,4,7,10)) +
  theme_bw() + 
  theme(text = element_text(size=FONT_SIZE), axis.text=element_text(size=FONT_SIZE),
        axis.line = element_line(),
        axis.text.y = element_text(size=FONT_SIZE, color = "gray30"),
        axis.title.y = element_text(size=FONT_SIZE, color = "black"),
        axis.text.x = element_text(size=FONT_SIZE, color="gray40"),
        legend.background = element_blank(), legend.position = "top",
        legend.text = element_text(size=FONT_SIZE-3))  +
  xlab("Schedule") + ylab("Number of Biopsies") + ggtitle("Progression speed: Slow")

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
