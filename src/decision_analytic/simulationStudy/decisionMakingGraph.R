meanOffset = by(scheduleResCombined, INDICES = scheduleResCombined$methodName, function(x){
  mean(x$offset[x$progression_time!=10])
})
meanNb = by(scheduleResCombined$nb, INDICES = scheduleResCombined$methodName, mean)
sdOffset = by(scheduleResCombined, INDICES = scheduleResCombined$methodName, function(x){
  sd(x$offset[x$progression_time!=10])
})
sdNb = by(scheduleResCombined$nb, INDICES = scheduleResCombined$methodName, sd)

medianOffset = by(scheduleResCombined, INDICES = scheduleResCombined$methodName, function(x){
  median(x$offset[x$progression_time!=10])
})
medianNb = by(scheduleResCombined$nb, INDICES = scheduleResCombined$methodName, median)

IQROffset = by(scheduleResCombined, INDICES = scheduleResCombined$methodName, function(x){
  IQR(x$offset[x$progression_time!=10])
})
IQRNb = by(scheduleResCombined$nb, INDICES = scheduleResCombined$methodName, IQR)

QRpt25Offset = by(scheduleResCombined, INDICES = scheduleResCombined$methodName, function(x){
  quantile(x$offset[x$progression_time!=10], probs = c(0.25))
})
QRpt25Nb = by(scheduleResCombined$nb, INDICES = scheduleResCombined$methodName, quantile, probs=0.25)

QRpt75Offset = by(scheduleResCombined, INDICES = scheduleResCombined$methodName, function(x){
  quantile(x$offset[x$progression_time!=10], probs = c(0.75))
})
QRpt75Nb = by(scheduleResCombined$nb, INDICES = scheduleResCombined$methodName, quantile, probs=0.75)

methodNames = levels(scheduleResCombined$methodName)

fixedScheduleIndices = c(1:3,14,15)
fixedScheduleLabels = c("Biopsy every\n1.5 years", "Annual\nbiopsies", 
                        "Biopsy every\n2 years", "Biopsy every\n3 years", "Biopsy every\n4 years")

FONT_SIZE=11
POINT_SIZE = 3
better_balance_intro = ggplot() + geom_label(aes(x=medianNb[fixedScheduleIndices], 
                                                   y=medianOffset[fixedScheduleIndices], 
                                                   label=fixedScheduleLabels),
                                               color='red2', size=2.75,
                                               nudge_x = c(0.6,-0.4,0.85,0.85,0.0), 
                                               nudge_y = c(0.175, 0.2, 0.125,0.125,0.125)) +
  geom_ribbon(aes(x=c(1:3,medianNb[fixedScheduleIndices]), ymin=0,
                  ymax=c(Inf,Inf,Inf,medianOffset[fixedScheduleIndices]), 
                  fill="Region of better balance in median number of biopsies and delay,\n than fixed/heuristic schedules"),alpha=0.2)+
  geom_point(aes(x=medianNb[fixedScheduleIndices], 
                 y=medianOffset[fixedScheduleIndices]), 
             color="red2", size=POINT_SIZE, shape=15) +
  geom_point(aes(x=1, y=0), color="dodgerblue4", size=POINT_SIZE) +
  geom_label(aes(x=1, y=0), label="Ideal\nschedule", 
             color="dodgerblue4", size=2.75, nudge_x = 0.65, nudge_y = 0.175) +
  geom_line(aes(x=medianNb[fixedScheduleIndices], 
                y=medianOffset[fixedScheduleIndices]), linetype="dashed", color="red2") +
  scale_fill_manual(name="", values="forestgreen") + 
  xlab("Median number of biopsies") + ylab("Median delay in detection of\ncancer progression (years)") +
  coord_cartesian(ylim=c(0,2)) +
  scale_x_continuous(breaks=1:10, limits = c(1,10)) +
  theme_bw() + 
  theme(text = element_text(size=FONT_SIZE), axis.text=element_text(size=FONT_SIZE),
        axis.line = element_line(), legend.background = element_blank(), legend.position = "bottom",
        legend.text = element_text(size=FONT_SIZE-3),
        plot.margin = margin(0, 0, 0, 0, "pt"))

ggsave(filename = "report/decision_analytic/mdm/latex/images/better_balance_intro.eps",
       plot=better_balance_intro, device=cairo_ps, height=5.5/1.333, width=5.5, dpi = 500)


## To make this plot using NPMLE plot
load("Rdata/decision_analytic/npmle_prias.Rdata")
survProb = 1 - cumsum(npmle_prias$pf)
survProb = c(1, survProb)

survIntervals = npmle_prias$intmap
survIntervals = cbind(c(0,0), survIntervals)

timePoints = as.numeric(survIntervals)[1:69]
survProbs = c(1,as.numeric(rep(survProb, each=2)))[1:69]

#Simulate failure times for 10000 patients
set.seed(2018)
survProbsPatients = runif(100000, min = 0, max = 1)
eventTimePatients = sapply(survProbsPatients, function(x){
  if(x < min(survProbs)){
    return(10)
  }
  
  minIndex = which.min(abs(x - survProbs))
  
  if(minIndex==1){
    return(runif(1, timePoints[1], timePoints[3]))
  }else{
    #print(survProbs[minIndex] == survProbs[minIndex + 1])
    return(runif(1, timePoints[minIndex], timePoints[minIndex + 1]))
  }
})
eventTimePatients[eventTimePatients>10] = 10

nbMedian_NPMLE = c(median(runFixedSchedule(eventTimePatients, seq(1.5, 10, 1.5))[,1]),
                   median(runFixedSchedule(eventTimePatients, seq(1, 10, 1))[,1]),
                   median(runFixedSchedule(eventTimePatients, seq(2, 10, 2))[,1]),
                   median(runFixedSchedule(eventTimePatients, seq(3, 10, 3))[,1]),
                   median(runFixedSchedule(eventTimePatients, seq(4, 10, 4))[,1]))

offsetMedian_NPMLE = c(median(runFixedSchedule(eventTimePatients, seq(1.5, 10, 1.5))[eventTimePatients!=10,2]),
                   median(runFixedSchedule(eventTimePatients, seq(1, 10, 1))[eventTimePatients!=10,2]),
                   median(runFixedSchedule(eventTimePatients, seq(2, 10, 2))[eventTimePatients!=10,2]),
                   median(runFixedSchedule(eventTimePatients, seq(3, 10, 3))[eventTimePatients!=10,2]),
                   median(runFixedSchedule(eventTimePatients, seq(4, 10, 4))[eventTimePatients!=10,2]))

fixedScheduleLabels = c("Biopsy every\n1.5 years", "Annual\nbiopsies", 
                        "Biopsy every\n2 years", "Biopsy every\n3 years", "Biopsy every\n4 years")

FONT_SIZE=11
POINT_SIZE = 3
better_balance_intro = ggplot() + geom_label(aes(x=nbMedian_NPMLE, 
                                                 y=offsetMedian_NPMLE, 
                                                 label=fixedScheduleLabels),
                                             color='red2', size=2.75,
                                             nudge_x = c(0.5,-0.4,0.7,0.7,0.0), 
                                             nudge_y = c(0.175, 0.2, 0.125,0.125,0.125)) +
  geom_ribbon(aes(x=c(1,2,nbMedian_NPMLE), ymin=0,
                  ymax=c(Inf,Inf,offsetMedian_NPMLE), 
                  fill="Region of better balance in median number of biopsies and delay,\n than fixed/heuristic schedules"),alpha=0.2)+
  geom_point(aes(x=nbMedian_NPMLE, 
                 y=offsetMedian_NPMLE), 
             color="red2", size=POINT_SIZE, shape=15) +
  geom_point(aes(x=1, y=0), color="dodgerblue4", size=POINT_SIZE) +
  geom_label(aes(x=1, y=0), label="Ideal\nschedule", 
             color="dodgerblue4", size=2.75, nudge_x = 0.5, nudge_y = 0.175) +
  geom_line(aes(x=nbMedian_NPMLE, 
                y=offsetMedian_NPMLE), linetype="dashed", color="red2") +
  scale_fill_manual(name="", values="forestgreen") + 
  xlab("Median number of biopsies") + ylab("Median delay in detection of\ncancer progression (years)") +
  coord_cartesian(ylim=c(0,2)) +
  scale_x_continuous(breaks=1:8, limits = c(1,8)) +
  theme_bw() + 
  theme(text = element_text(size=FONT_SIZE), axis.text=element_text(size=FONT_SIZE),
        axis.line = element_line(), legend.background = element_blank(), legend.position = "bottom",
        legend.text = element_text(size=FONT_SIZE-3),
        plot.margin = margin(0, 0, 0, 0, "pt"))

ggsave(filename = "report/decision_analytic/mdm/latex/images/better_balance_intro.eps",
       plot=better_balance_intro, device=cairo_ps, height=5.5/1.333, width=5.5, dpi = 500)


########################
fixedScheduleIndices = c(1:3,14,15)
fixedScheduleLabels = c("Biopsy every\n1.5 years", "Annual\nbiopsies", 
                        "Biopsy every\n2 years", "Biopsy every\n3 years", "Biopsy every\n4 years")

priasIndex = 7
priasLabel = "PRIAS"

persScheduleIndices = c(8,9,11,13)
persScheduleLabels = c("Risk: 10%", "Risk: 15%", "Risk: 5%", "Risk: F1") 

FONT_SIZE=11
POINT_SIZE = 3

better_balance_results= ggplot() + geom_label(aes(x=medianNb[fixedScheduleIndices], 
                                                 y=medianOffset[fixedScheduleIndices], 
                                                 label=fixedScheduleLabels),
                                             color='red2', size=2.75,
                                             nudge_x = c(0.6,-0.4,0.85,0.85,0.0), 
                                             nudge_y = c(0.175, 0.2, 0.125,0.125,0.125)) +
  geom_ribbon(aes(x=c(1:3,medianNb[fixedScheduleIndices]), ymin=0,
                  ymax=c(Inf,Inf,Inf,medianOffset[fixedScheduleIndices]), 
                  fill="Region of better balance in median\n number of biopsies and delay,\n than fixed/heuristic schedules"),alpha=0.2)+
  geom_line(aes(x=medianNb[fixedScheduleIndices], 
                y=medianOffset[fixedScheduleIndices]), linetype="dashed", color="red2") +
  
  geom_point(aes(x=medianNb[fixedScheduleIndices], 
                 y=medianOffset[fixedScheduleIndices], shape="Fixed", color="Fixed"),
             size=POINT_SIZE) +
  
  geom_point(aes(x=medianNb[persScheduleIndices], 
                 y=medianOffset[persScheduleIndices], shape="Personalized", color="Personalized"), 
            size=POINT_SIZE) +
  geom_label(aes(x=medianNb[persScheduleIndices], 
                 y=medianOffset[persScheduleIndices], 
                 label=persScheduleLabels),
             color='forestgreen', size=2.75,
             nudge_y = c(-0.15, 0.15, -0.15, -0.15)) +
  geom_line(aes(x=medianNb[persScheduleIndices], 
                y=medianOffset[persScheduleIndices]), linetype="dashed", color="forestgreen") +
  
  geom_point(aes(x=medianNb[priasIndex], y=medianOffset[priasIndex], shape="PRIAS", color="PRIAS"), 
             size=POINT_SIZE, alpha=0.8) +
  geom_label(aes(x=medianNb[priasIndex],
                 y=medianOffset[priasIndex]), label=priasLabel, 
             color="dodgerblue4", size=2.75, nudge_x = 0.3, nudge_y = 0.15) +
  
  scale_shape_manual(name="", values=c(15,17,18)) + 
  scale_color_manual(name="", values=c("red2", "forestgreen", "dodgerblue4")) + 
  
  scale_fill_manual(name="", values="forestgreen") + 
  xlab("Median number of biopsies") + ylab("Median delay in detection of\ncancer progression (years)") +
  coord_cartesian(ylim=c(0,2)) +
  scale_x_continuous(breaks=1:10, limits = c(1,10)) +
  theme_bw() + 
  theme(text = element_text(size=FONT_SIZE), axis.text=element_text(size=FONT_SIZE),
        axis.line = element_line(), legend.background = element_blank(), legend.position = "bottom",
        legend.text = element_text(size=FONT_SIZE-3),
        plot.margin = margin(0, 0, 0, 0, "pt"))

ggsave(filename = "report/decision_analytic/mdm/latex/images/better_balance_results.eps",
       plot=better_balance_results, device=cairo_ps, height=5.5/1.333, width=5.5, dpi = 500)


###### 
# Comparing PRIAS and 10%
######
scheduleResCombined$yearProgression = ceiling(scheduleResCombined$progression_time)
scheduleResCombined$yearProgression = sapply(scheduleResCombined$progression_time, function(t){
  if(t<=3.5){
    return("Fast")
  }else if(t==10){
    "Slow"
  }else{
    "Intermediate"
  }
})

by(scheduleResCombined, INDICES = scheduleResCombined$yearProgression, function(x){
  by(x$nb, INDICES = x$methodName, median)
})

ggplot(data=scheduleResCombined) + geom_point(aes(x=methodName, y=nb, color=yearProgression, group=yearProgression))
