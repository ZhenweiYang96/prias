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
POINT_SIZE = 2
better_balance_intro = ggplot() + geom_label(aes(x=medianNb[fixedScheduleIndices], 
                                                   y=medianOffset[fixedScheduleIndices], 
                                                   label=fixedScheduleLabels),
                                               color='red2', size=2.75,
                                               nudge_x = c(0.6,-0.4,0.85,0.85,0.0), 
                                               nudge_y = c(0.175, 0.2, 0.125,0.125,0.125)) +
  geom_ribbon(aes(x=c(1:3,medianNb[fixedScheduleIndices]), ymin=0,
                  ymax=c(Inf,Inf,Inf,medianOffset[fixedScheduleIndices]), 
                  fill="Region of better balance in median number of biopsies and delay than fixed schedules"),alpha=0.2)+
  geom_point(aes(x=medianNb[fixedScheduleIndices], 
                 y=medianOffset[fixedScheduleIndices]), color="red2", size=POINT_SIZE) +
  geom_point(aes(x=1, y=0), color="dodgerblue4", size=POINT_SIZE) +
  geom_label(aes(x=1, y=0), label="Ideal\nschedule", 
             color="dodgerblue4", size=2.75, nudge_x = 0.65, nudge_y = 0.175) +
  geom_line(aes(x=medianNb[fixedScheduleIndices], 
                y=medianOffset[fixedScheduleIndices]), linetype="dashed", color="red2") +
  scale_fill_manual(name="", values="forestgreen") + 
  xlab("Median Number of biopsies") + ylab("Median delay in detection of\ncancer progression (years)") +
  coord_cartesian(ylim=c(0,2)) +
  scale_x_continuous(breaks=1:10, limits = c(1,10)) +
  theme_bw() + 
  theme(text = element_text(size=FONT_SIZE), axis.text=element_text(size=FONT_SIZE),
        axis.line = element_line(), legend.background = element_blank(), legend.position = "bottom",
        legend.text = element_text(size=FONT_SIZE-3),
        plot.margin = margin(0, 0, 0, 0, "pt"))

ggsave(filename = "report/decision_analytic/mdm/latex/images/better_balance_intro.eps",
       plot=better_balance_intro, device=cairo_ps, height=5.5/1.333, width=5.5, dpi = 500)



mean_nb_offset = ggplot() + geom_label(aes(x=meanNb, y=meanOffset, label=methodNames)) +
  xlab("Mean Number of biopsies") + ylab("Mean delay in detection of cancer progression (years)") +
  ylim(0, 2) +
  scale_x_continuous(breaks=1:7, limits = c(1,7)) +
  theme_bw() + 
  theme(text = element_text(size=FONT_SIZE), axis.text=element_text(size=FONT_SIZE),
        axis.line = element_line())

ggsave(filename = "report/decision_analytic/mdm/latex/images/mean_nb_offset.eps",
       plot=mean_nb_offset, device=cairo_ps, height=5.5/1.333, width=5.5, dpi = 500)

#This plot the advantage of personalized over fixed schedules
fixedIndices = c(1,2,3,16)
persIndices = c(8,9,10,11)

fixedIndices=7
persIndices=8

ggplot() + 
  #geom_line(aes(x=medianNb[fixedIndices], y=medianOffset[fixedIndices], color='Fixed'), linetype="dotted") +
  geom_point(aes(x=medianNb[fixedIndices], y=medianOffset[fixedIndices], 
                 label=methodNames[fixedIndices], color='Fixed')) +
  geom_errorbarh(aes(y=medianOffset[fixedIndices], xmin=QRpt25Nb[fixedIndices], 
                    xmax=QRpt75Nb[fixedIndices], color='Fixed'), width=0.1) +
  geom_line(aes(x=medianNb[persIndices], y=medianOffset[persIndices], color='Personalized'), linetype="dotted") +
  geom_point(aes(x=medianNb[persIndices], y=medianOffset[persIndices],
                 label=methodNames[persIndices], color='Personalized')) +
  geom_errorbarh(aes(y=medianOffset[persIndices], xmin=QRpt25Nb[persIndices], 
                    xmax=QRpt75Nb[persIndices], color='Personalized'), width=0.1) +
  xlab("Median Number of biopsies") + ylab("Median delay in detection of cancer progression (years)") +
  ylim(0, 2.5) +
  scale_x_continuous(breaks=1:11, limits = c(1,11)) +
  theme_bw() + 
  theme(text = element_text(size=FONT_SIZE), axis.text=element_text(size=FONT_SIZE),
        axis.line = element_line()) + coord_flip()

ggplot() + geom_line(aes(x=medianNb[fixedIndices], y=medianOffset[fixedIndices], color='Fixed'), linetype="dotted") +
  geom_point(aes(x=medianNb[fixedIndices], y=medianOffset[fixedIndices], 
                 label=methodNames[fixedIndices], color='Fixed')) +
  geom_errorbar(aes(x=medianNb[fixedIndices], ymin=QRpt25Offset[fixedIndices], 
                    ymax=QRpt75Offset[fixedIndices], color='Fixed'), width=0.25) +
  geom_line(aes(x=medianNb[persIndices], y=medianOffset[persIndices], color='Personalized'), linetype="dotted") +
  geom_point(aes(x=medianNb[persIndices], y=medianOffset[persIndices],
                 label=methodNames[persIndices], color='Personalized')) +
  geom_errorbar(aes(x=medianNb[persIndices], ymin=QRpt25Offset[persIndices], 
                    ymax=QRpt75Offset[persIndices], color='Personalized'), width=0.25) +
  xlab("Median Number of biopsies") + ylab("Median delay in detection of cancer progression (years)") +
  ylim(0, 2.5) +
  scale_x_continuous(breaks=1:11, limits = c(1,11)) +
  theme_bw() + 
  theme(text = element_text(size=FONT_SIZE), axis.text=element_text(size=FONT_SIZE),
        axis.line = element_line())

#Make a plot with year of progression on X axis, color for schedule, mean delay on Y axis, 
#and size for number of biopsies
scheduleResCombined$yearProgression = ceiling(scheduleResCombined$progression_time)
scheduleResCombined$yearProgression = sapply(scheduleResCombined$progression_time, function(x){
  if(x<=3.5){
    "Fast"
  }else if(x==10){
    "Slow"
  }else{
    "Intermediate"
  }
})

yearlyRes = by(scheduleResCombined, INDICES = scheduleResCombined$yearProgression, FUN = function(yearResCombined){
  res = by(yearResCombined, INDICES = yearResCombined$methodName, FUN = function(methodRes){
    c(median(methodRes$nb), median(methodRes$offset[methodRes$progression_time!=10]), nrow(methodRes))
  })
  
  resMatrix = do.call(rbind,res)
  resDf= data.frame(resMatrix)
  colnames(resDf) = c("nb", "offset", "nPatients")
  resDf$methodName = rownames(resMatrix)
  resDf$yearProgression  = yearResCombined$yearProgression[1]
  return(resDf)
})

yearlyResDf = do.call(rbind, yearlyRes)
yearlyResDf$offset[is.na(yearlyResDf$offset)] = runif(n = sum(is.na(yearlyResDf$offset)),-1, 0)
ggplot(data=yearlyResDf[yearlyResDf$methodName %in% c("Annual", "PRIAS", "Risk (10%)"),]) + 
  geom_point(aes(x=yearProgression, y=offset, color=methodName, size=nb)) +
  #scale_x_continuous(breaks = 1:10, limits = c(1,10)) + theme_bw() + 
  theme(text = element_text(size=FONT_SIZE), axis.text=element_text(size=FONT_SIZE),
        axis.line = element_line()) + ylab("Delay (years)")

