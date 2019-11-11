library(ggplot2)
library(ggpubr)

MAX_FOLLOW_UP = 10

fixed_df_10 = fixed_df[fixed_df$MAX_FOLLOW_UP==MAX_FOLLOW_UP,]
fixed_df_10$ScheduleType="Fixed"
fixed_df_10$ScheduleType[fixed_df_10$schedule=="PRIAS"] = "PRIAS"

pers_10 = personalized_df[personalized_df$MAX_FOLLOW_UP==MAX_FOLLOW_UP,]
pers_10[,2] = paste(pers_10[,2], pers_10[,3], sep = "-")
pers_10 = pers_10[,-3]
pers_10$ScheduleType="Personalized"

colnames(pers_10)[2] = 'schedule'

plotDf = rbind(fixed_df_10, pers_10)
plotDf = plotDf[!(plotDf$schedule %in% c("Triennial", "Triennial_1")),]
plotDf$schedule = factor(plotDf$schedule)

#levels(plotDf$schedule)[5:7] = paste("Max delay", c(0.5, 1, Inf), "year")
levels(plotDf$schedule)[5:9] = paste0("Max threshold ", c(0.05, 0.1, 0.15, 0.25, 0.5)*100, "%")

a = ggplot(data=plotDf[plotDf$progression_time == plotDf$MAX_FOLLOW_UP,]) + 
  geom_boxplot(aes(x=schedule, y=nb, color=ScheduleType), outlier.shape = NA) + 
  ylab("#Biopsies") +
  scale_y_continuous(breaks=0:MAX_FOLLOW_UP, limits = c(0, MAX_FOLLOW_UP))+
  theme_bw() + theme(text = element_text(size=14), axis.text.x = element_text(angle = 90, hjust = 1)) + 
  ggtitle("Non Progressors")

b = ggplot(data=plotDf[plotDf$progression_time < plotDf$MAX_FOLLOW_UP,]) + 
  geom_boxplot(aes(x=schedule, y=delay, color=ScheduleType), outlier.shape = NA) + 
  ylab("Delay (years)") +
  scale_y_continuous(breaks=0:5, limits = c(0, 5))+
  theme_bw() + theme(text = element_text(size=14), axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Progressors")

c = ggplot(data=plotDf[plotDf$progression_time < plotDf$MAX_FOLLOW_UP,]) + 
  geom_boxplot(aes(x=schedule, y=nb, color=ScheduleType), outlier.shape = NA) + 
  ylab("#Biopsies") +
  scale_y_continuous(breaks=0:MAX_FOLLOW_UP, limits = c(0, MAX_FOLLOW_UP))+
  theme_bw() + theme(text = element_text(size=14), axis.text.x = element_text(angle = 90, hjust = 1))+ 
  ggtitle("Progressors")

ggpubr::ggarrange(a,c,b, ncol = 3, nrow = 1, common.legend = T)
