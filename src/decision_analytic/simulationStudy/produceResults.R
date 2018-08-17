#First load all the results and extract the actual schedule results
library(doParallel)

resFiles = list.files("/home/a_tomer/Results/final_res_2nd_paper/", full.names = T)

prias_real_bokhorst = foreach(i=2:length(resFiles), .combine="rbind") %do%{
  load(resFiles[i])
  
  print(resFiles[i])
  
  seed = jointModelData$seed
  set.seed(seed)
  
  real_f1 = jointModelData$testData$scheduleResults[jointModelData$testData$scheduleResults$methodName==KAPPAF1Score,]
  real_f1$methodName = "Risk (F1) Real"
  
  real_f1[,c("nb", "offset")] = runDynRiskGRSchedule(jointModelData, "F1score")
  jointModelData$testData$scheduleResults = rbind(jointModelData$testData$scheduleResults, real_f1)
  
  save(jointModelData, file=resFiles[i])
  ret = jointModelData$testData$scheduleResults
  rm(jointModelData)
  return(ret)
}

scheduleResCombined = prias_real_bokhorst[prias_real_bokhorst$methodName %in% c(ANNUAL, PRIAS, KAPPApt85, KAPPApt95, KAPPAF1Score, "Risk (F1) Real"),]
scheduleResCombined$methodName = droplevels(scheduleResCombined$methodName)

scheduleResCombined$nb[scheduleResCombined$offset < 0] = scheduleResCombined$nb[scheduleResCombined$offset<0] + 1
scheduleResCombined$offset[scheduleResCombined$offset < 0] = 10 - scheduleResCombined$progression_time[scheduleResCombined$offset<0]

levels(scheduleResCombined$methodName)[5] = "Risk (F1): Dt=0.5"
levels(scheduleResCombined$methodName)[6] = "Risk (F1): Dt=s-t"

ggplot(data=scheduleResCombined[scheduleResCombined$progression_time==10,]) + 
  geom_boxplot(aes(x=methodName, y=nb), outlier.shape = NA) + coord_flip(ylim = c(0,10)) +
  theme(text = element_text(size=17.5)) + xlab("Schedule") + ylab("Number of Biopsies") +
  geom_hline(yintercept = 3, linetype="dashed", color='red') + 
  geom_hline(yintercept = 0, linetype="dashed", color='black')

ggsave(.Last.value, filename = "biopsycens.eps", width = 8.27)

gfastNb = ggplot(data=scheduleResCombined[scheduleResCombined$progression_time>0 & scheduleResCombined$progression_time<=3.5,]) + 
  geom_boxplot(aes(x=methodName, y=nb), outlier.shape = NA) + coord_flip(ylim = c(0, 3)) +
  theme(text = element_text(size=17.5)) + xlab("Schedule") + ylab("Number of Biopsies") +
  geom_hline(yintercept = 3, linetype="dashed", color='red')

gfastOffset = ggplot(data=scheduleResCombined[scheduleResCombined$progression_time>0 & scheduleResCombined$progression_time<=3.5,]) + 
  geom_boxplot(aes(x=methodName, y=offset), outlier.shape = NA) + coord_flip(ylim = c(0,4.5)) +
  theme(text = element_text(size=17.5)) + xlab("Schedule") + ylab("Undetected time in years (last biopsy time - true GR time)") +
  geom_hline(yintercept = 2, linetype="dashed", color='red') + geom_hline(yintercept = 0, linetype="dashed", color='black')

ggpubr::ggarrange(gfastNb, gfastOffset, nrow=2, labels = "AUTO")
ggsave(.Last.value, file="biopsy3pt5.eps", height = 9, width = 8.27)

gslowNb = ggplot(data=scheduleResCombined[scheduleResCombined$progression_time>3.5 & scheduleResCombined$progression_time<10,]) + 
  geom_boxplot(aes(x=methodName, y=nb), outlier.shape = NA) + coord_flip(ylim = c(0, 10)) +
  theme(text = element_text(size=17.5)) + xlab("Schedule") + ylab("Number of Biopsies") +
  geom_hline(yintercept = 3, linetype="dashed", color='red')

gslowOffset = ggplot(data=scheduleResCombined[scheduleResCombined$progression_time>3.5 & scheduleResCombined$progression_time<10,]) + 
  geom_boxplot(aes(x=methodName, y=offset), outlier.shape = NA) + coord_flip(ylim = c(0,6)) +
  theme(text = element_text(size=17.5), axis.title = element_text(size=17.5)) + xlab("Schedule") + ylab("Undetected time in years (last biopsy time - true GR time)") +
  geom_hline(yintercept = 2, linetype="dashed", color='red') + geom_hline(yintercept = 0, linetype="dashed", color='black')

ggpubr::ggarrange(gslowNb, gslowOffset, nrow=2, labels = "AUTO")
ggsave(ggpubr::ggarrange(gslowNb, gslowOffset, nrow=2, labels = "AUTO"), 
       file="biopsymedium.eps", height = 9, width = 8.27)


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
