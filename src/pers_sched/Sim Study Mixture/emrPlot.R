install.packages("ggpubr")
library(ggpubr)

load("Rdata/Gleason as event/Sim Study/biopsyResultsAll.Rdata")
biopsyResultsAll = biopsyResultsAll[biopsyResultsAll$methodName %in% c("Exp. GR Time", "Med. GR Time", "Dyn. Risk GR", "Annual", "Hybrid", "PRIAS"),]
biopsyResultsAll$methodName = droplevels(biopsyResultsAll$methodName)
levels(biopsyResultsAll$methodName)
levels(biopsyResultsAll$methodName) = c("Exp. Time", "Dyn. risk", "Med. Time", "PRIAS", "Annual", "Hybrid")
biopsyResultsAll$Schedule.Type = ifelse(biopsyResultsAll$methodName %in% c("PRIAS", "Annual"), "Fixed", "Personalized")

source("../JMBayes/Anirudh/dev/multiplot.R")

gnb_all = ggplot(data = biopsyResultsAll) +
  geom_boxplot(aes(reorder(methodName, nb, FUN=median), nb, color=Schedule.Type), outlier.shape = NA) +
  coord_flip(ylim = c(1, 14)) + scale_y_continuous(breaks=seq(1, 14, 3)) +
  theme(text = element_text(size=14), axis.text=element_text(size=14), legend.text = element_text(size=14)) +
  ylab("Number of biopsies") + xlab("Schedule") + ggtitle("All patients")

goffset_all = ggplot(data = biopsyResultsAll) +
  geom_boxplot(aes(reorder(methodName, nb, FUN=median), offset*12, color=Schedule.Type), outlier.shape = NA) +
  theme(text = element_text(size=14), axis.text=element_text(size=14), legend.text = element_text(size=14)) +
  coord_flip(ylim=c(1,45)) + 
  ylab("Delay in detection of progression (months)") + xlab("Schedule") + ggtitle("All patients")

gnb_slow = ggplot(data = biopsyResultsAll[biopsyResultsAll$weibullScale==6,]) +
  geom_boxplot(aes(reorder(methodName, nb, FUN=median), nb, color=Schedule.Type), outlier.shape = NA) +
  theme(text = element_text(size=14), axis.text=element_text(size=14), legend.text = element_text(size=14)) +
  coord_flip(ylim = c(1, 14)) + scale_y_continuous(breaks=seq(1, 14, 3)) +
  ylab("Number of biopsies") + xlab("Schedule") + ggtitle("Slowly progressing patients")

goffset_slow = ggplot(data = biopsyResultsAll[biopsyResultsAll$weibullScale==6,]) +
  theme(text = element_text(size=14), axis.text=element_text(size=14), legend.text = element_text(size=14)) +
  geom_boxplot(aes(reorder(methodName, nb, FUN=median), offset*12, color=Schedule.Type), outlier.shape = NA) +
  coord_flip(ylim=c(1,45)) + 
  ylab("Delay in detection of progression (months)") + xlab("Schedule")  + ggtitle("Slowly progressing patients")

meanNb_all = c(by(data = biopsyResultsAll$nb, INDICES = biopsyResultsAll$methodName, mean))
meanOffset_all = c(by(data = biopsyResultsAll$offset*12, INDICES = biopsyResultsAll$methodName, mean))

meanNb_slow = c(by(data = biopsyResultsAll$nb[biopsyResultsAll$weibullScale==6], INDICES = biopsyResultsAll$methodName[biopsyResultsAll$weibullScale==6], mean))
meanOffset_slow = c(by(data = biopsyResultsAll$offset[biopsyResultsAll$weibullScale==6]*12, INDICES = biopsyResultsAll$methodName[biopsyResultsAll$weibullScale==6], mean))


gnbvsoffset_all = qplot(x = meanNb_all, 
                        y = meanOffset_all, 
                        label=levels(biopsyResultsAll$methodName),
                        xlab="Mean number of biopsies", ylab="Mean biopsy offset (months)")+ 
  coord_cartesian(xlim=c(0,5.5), ylim=c(0, 16)) + 
  geom_point() +  geom_label(size=3.5, nudge_y = -0.6)

gnbvsoffset_slow = qplot(x = meanNb_slow, 
                        y = meanOffset_slow, 
                        label=levels(biopsyResultsAll$methodName),
                        xlab="Mean number of biopsies", ylab="Mean biopsy offset (months)")+ 
  coord_cartesian(xlim=c(0,6.5), ylim=c(0, 11)) + 
  geom_point() +  geom_label(size=3.5, nudge_y = -0.6)


resPlots = ggarrange(gnb_all, goffset_all, gnb_slow, goffset_slow,
                     labels=c("A", "B", "C","D"),
          ncol = 2, nrow = 2, common.legend = T, legend = "bottom")

ggsave(file="../../docs/Presentation/Conference 2018/EMR Extended Abstract/resPlots.eps", plot = resPlots,width=10, height = 8.27)
