training.prias.id = prias.id[!(prias.id$P_ID %in% c(3174, 2340, 911)),]

###########################################################################
survModel = survreg(Surv(progression_time_start, progression_time_end, type = "interval2") ~ 
                       I(Age - 70) +  I((Age - 70)^2), data = prias.id, model = TRUE)
survModel.training= survreg(Surv(progression_time_start, progression_time_end, type = "interval2") ~ 
                              I(Age - 70) +  I((Age - 70)^2), data = training.prias.id, model = TRUE)

#interval censoring NPMLE. Takes time to run. beware
npmle_prias = icfit(Surv(progression_time_start,progression_time_end,type="interval2")~1, data=prias.id)
#The density function is given by npmle_prias$pf
survProb = 1 - cumsum(npmle_prias$pf)
survProb = c(1, survProb)

survIntervals = npmle_prias$intmap
survIntervals = cbind(c(0,0), survIntervals)

timePoints = as.numeric(survIntervals)
survProbs = c(1,as.numeric(rep(survProb, each=2)))[1:length(timePoints)]

FONT_SIZE=11
npmle_plot = ggplot() + geom_line(aes(x=timePoints, y=1-survProbs)) +  
  coord_cartesian(xlim=c(0,10)) + 
  theme_bw() +
  theme(text = element_text(size=FONT_SIZE), axis.text=element_text(size=FONT_SIZE),
        axis.line = element_line(), 
        plot.margin = margin(0, 0, 0, 0, "pt")) + 
  scale_y_continuous(breaks = seq(0, 1, 0.25), labels = paste0(seq(0, 1, 0.25)*100, "%"),
                    limits = c(0,1)) + 
  ylab("Cumulative risk of cancer progression (%)") +
  xlab("Follow-up time (years)")

ggsave(filename = "report/decision_analytic/mdm/latex/images/npmle_plot.eps",
       plot=npmle_plot, device=cairo_ps, height=5.5/1.333, width=5.5, dpi = 500)

#############################################################
# By assuming right censoring
#############################################################
prias.id.rightCens = prias.id
prias.id.rightCens$progressed = ifelse(prias.id$progressed>0, 1, 0)
prias.id.rightCens$progression_time = ifelse(prias.id$progression_time_end==Inf, 
                                             prias.id$progression_time_start, 
                                             0.5 * (prias.id$progression_time_start+prias.id$progression_time_end))

kmfit = survfit(Surv(progression_time, progressed)~1, conf.type="log-log", data=prias.id.rightCens)
survminer::ggsurvplot(kmfit, risk.table = T,break.time.by = 1, xlab = "Time (years)", ylim = c(0.5,1), conf.int = T)


survModel_rightCens = coxph(Surv(progression_time, progressed) ~ I(Age - 70) +  I((Age - 70)^2), 
                                    data=prias.id.rightCens, x = T, model = T)
#

training.prias.id.rightCens = training.prias.id
training.prias.id.rightCens$progressed = ifelse(training.prias.id$progressed>0, 1, 0)
training.prias.id.rightCens$progression_time = ifelse(training.prias.id$progression_time_end==Inf, 
                                                      training.prias.id$progression_time_start, 
                                                      training.prias.id$progression_time_end)

survModel.training_rightCens = coxph(Surv(progression_time, progressed) ~ I(Age - 70) +  I((Age - 70)^2), 
                            data=training.prias.id.rightCens, x = T, model = T)
