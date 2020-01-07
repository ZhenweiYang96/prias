library(JMbayes)
library(splines)
library(ggplot2)
library(ggpubr)

load("Rdata/lastpaper/fitted_model/mvJoint_dre_psa_2knots_quad_age_light.Rdata")
load("Rdata/gap3/PRIAS_2019/cleandata.Rdata")
source("src/lastpaper/prediction_psa_dre.R")
source("src/lastpaper/scheduleCreator.R")

SUCCESS_COLOR = 'forestgreen'
DANGER_COLOR = 'red'
WARNING_COLOR = 'darkorange'
THEME_COLOR = 'dodgerblue4'
MAX_FOLLOW_UP = 10
POINT_SIZE = 1.75
FONT_SIZE = 11
LABEL_SIZE = 1.25

set.seed(2019)
patient_df = prias_long_final[prias_long_final$P_ID == 101 & 
                                prias_long_final$year_visit<=5,]
patient_df$palpable_dre = as.numeric(patient_df$palpable_dre)
patient_df$palpable_dre[1] = 1
patient_df$year_visit[nrow(patient_df)] = 5
current_visit_time = max(patient_df$year_visit)

last_biopsy_time = 3.5

surv_prob_times = seq(last_biopsy_time, MAX_FOLLOW_UP, length.out = 100)
future = getExpectedFutureOutcomes(object = mvJoint_dre_psa_2knots_quad_age_light, 
                                   patient_data = patient_df, 
                                   latest_survival_time = last_biopsy_time, 
                                   long_predict_times = patient_df$year_visit,
                                   survival_predict_times = surv_prob_times)

risk_probs = 1 - future$predicted_surv_prob

patient_df$truePSA = rowMeans(future$predicted_psa)
patient_df$trueProbPalpableDRE = rowMeans(future$predicted_dre_prob)

min_x = -0.4

DRE_PSA_Y_GAP = 0.1
minYLeft = 0
maxYleft = 2 + DRE_PSA_Y_GAP * 2

dreDs = patient_df[!is.na(patient_df$palpable_dre), c("year_visit", "palpable_dre", "trueProbPalpableDRE")]
psaDs = patient_df[!is.na(patient_df$psa), c("year_visit", "log2psaplus1", "truePSA")]
maxPSA = max(psaDs[,-1], na.rm = T)
minPSA = min(psaDs[,-1], na.rm = T)

newMinPSA = (maxYleft/2 + DRE_PSA_Y_GAP)
newMaxPSA = maxYleft

psaDs[,-1] = ((psaDs[,-1] - minPSA) * (newMaxPSA - newMinPSA))/ (maxPSA - minPSA) + newMinPSA

xTicks = c(0, last_biopsy_time, current_visit_time)
xLabels = c("0\n(Start AS)", paste0("t = ", round(last_biopsy_time,1), "\n(Latest\nbiopsy)"),
            paste0("s = ", round(current_visit_time,1), "\n(Current\nvisit)"))

drebreaks = seq(0, maxYleft/2 - DRE_PSA_Y_GAP, length.out = 3)
drelabels = paste0(drebreaks * 100, "%")
psabreaks = seq(maxYleft/2 + DRE_PSA_Y_GAP, maxYleft, length.out = 3)
psalabels = round(seq(minPSA, maxPSA, length.out = 3),1)

mean_risk_probs = c(0, apply(risk_probs, 1, mean))
upper_risk_probs = c(0, apply(risk_probs, 1, quantile, probs=0.975))
lower_risk_probs = c(0,apply(risk_probs, 1, quantile, probs=0.025))

scaled_mean_risk_probs = mean_risk_probs * (maxYleft - minYLeft) + minYLeft
scaled_upper_risk_probs = upper_risk_probs * (maxYleft - minYLeft) + minYLeft
scaled_lower_risk_probs = lower_risk_probs * (maxYleft - minYLeft) + minYLeft

riskAxisBreaks = c(0, 0.25, 0.5, 0.75, 1) * (maxYleft - minYLeft) + minYLeft
riskAxisLabels = c("0%", "25%", "50%", "75%", "100%")

baseggplot = ggplot() + 
  theme_bw() + 
  theme(plot.margin = margin(0,0,0,0, "mm"), 
        axis.text = element_text(size=FONT_SIZE),
        axis.title = element_text(size=FONT_SIZE),
        plot.title = element_text(size=FONT_SIZE),
        axis.text.y.left = element_text(size=FONT_SIZE, color=THEME_COLOR),
        axis.title.y.left = element_text(size=FONT_SIZE, color=THEME_COLOR),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size=FONT_SIZE-2))

baseggplot_no_xticks = baseggplot + 
  theme(axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())

risk_plot = baseggplot_no_xticks +
  geom_line(aes(x=surv_prob_times, y=scaled_mean_risk_probs), color=DANGER_COLOR) +
  geom_line(aes(x=surv_prob_times, y=scaled_mean_risk_probs), color=DANGER_COLOR) +
  geom_ribbon(aes(x=surv_prob_times, ymin=scaled_lower_risk_probs,
                  ymax=scaled_upper_risk_probs), alpha=0.15, fill=DANGER_COLOR) + 
  geom_segment(aes(x=-Inf, xend=current_visit_time, 
                   y=maxYleft/2, yend=maxYleft/2), linetype="solid") +
  geom_vline(xintercept = last_biopsy_time, color=SUCCESS_COLOR) + 
  geom_vline(xintercept = current_visit_time, color='black', linetype='dashed') + 
  geom_point(data = psaDs, size=POINT_SIZE, color=THEME_COLOR,
             aes(x = year_visit, y=log2psaplus1, shape="Observed PSA")) +
  geom_line(data = psaDs,  color=THEME_COLOR,
            aes(x = year_visit, y=truePSA, linetype="Fitted PSA")) +
  geom_point(data = dreDs, size=POINT_SIZE, color=THEME_COLOR,
             aes(x = year_visit, y=palpable_dre, shape="Observed DRE")) +
  geom_line(data = dreDs, color=THEME_COLOR,
            aes(x = year_visit, y=trueProbPalpableDRE, linetype="Fitted DRE")) +
  scale_shape_manual(name="",
                     labels= c(expression(atop('Observed', 'DRE')), expression(atop('Observed', 'log'[2]*'(PSA + 1)'))),
                     values = c(17, 16)) +
  scale_linetype_manual(name="",
                        labels= c(expression(atop('Fitted', 'Pr (Palpable DRE)')), expression(atop('Fitted', 'log'[2]*'(PSA + 1)'))),
                        values = c("dotted", "dashed")) +
  ylab(expression('Pr (DRE Palpable)   '*'log'[2]*'(PSA + 1)')) +
  #ylab('DRE              PSA') +
  xlab("Follow-up time (years)") + 
  xlim(min_x,MAX_FOLLOW_UP) +
  theme(axis.text.y.right = element_text(size=FONT_SIZE, color=DANGER_COLOR),
        axis.title.y.right = element_text(size=FONT_SIZE, color=DANGER_COLOR)) +
  scale_y_continuous(limits = c(minYLeft, maxYleft),
                     breaks = c(drebreaks, psabreaks),
                     labels = c(drelabels, psalabels),
                     sec.axis = sec_axis(trans=~., 
                                         breaks= riskAxisBreaks,
                                         labels = riskAxisLabels,
                                         name = "Cumulative-risk of progression")) 

print(risk_plot)

prias_schedule = getPRIASSchedule(object = mvJoint_dre_psa_2knots_quad_age_light,
                                  patient_data = patient_df, 
                                  cur_visit_time = current_visit_time, 
                                  last_biopsy_time = last_biopsy_time,
                                  M=500, horizon = MAX_FOLLOW_UP)

risk_10_perc_schedule = getFixedRiskBasedScheduleNoGrid(object = mvJoint_dre_psa_2knots_quad_age_light,
                                    patient_data = patient_df,
                                    cur_visit_time = current_visit_time,
                                    last_biopsy_time = last_biopsy_time,
                                    risk_threshold = 0.1, min_biopsy_gap = 1,
                                    M = 500, horizon = MAX_FOLLOW_UP)

t1 = Sys.time()
risk_automatic_schedule = getAutomaticRiskBasedScheduleNoGrid(object = mvJoint_dre_psa_2knots_quad_age_light,
                                                        patient_data = patient_df,
                                                        cur_visit_time = current_visit_time,
                                                        last_biopsy_time = last_biopsy_time,
                                                        min_biopsy_gap = 1,
                                                        M = 500, horizon = MAX_FOLLOW_UP)
t2 = Sys.time()
print(t2-t1)


schedule_df = data.frame(name=character(), number=numeric(),times=vector())
schedule_df = rbind(schedule_df, data.frame(name="Annual", 
                                            number=1,
                                            times=seq(current_visit_time, MAX_FOLLOW_UP, by = 1)))
schedule_df = rbind(schedule_df, data.frame(name="PRIAS", 
                                            number=2,
                                            times=prias_schedule))
schedule_df = rbind(schedule_df, data.frame(name="Risk: Auto", 
                                            number=3,
                                            times=risk_automatic_schedule$practical_biopsy_times))
schedule_df = rbind(schedule_df, data.frame(name="Risk: 10%", 
                                            number=4,
                                            times=risk_10_perc_schedule))

schedule_plot = baseggplot +
  geom_vline(xintercept = last_biopsy_time, color=SUCCESS_COLOR) + 
  geom_vline(xintercept = current_visit_time, color='black', linetype='dashed') + 
  geom_segment(aes(x=rep(current_visit_time,4), xend=rep(MAX_FOLLOW_UP,4), y=1:4, yend=1:4), color=SUCCESS_COLOR, linetype='dotted')+
  geom_label(data=schedule_df,
             aes(x=times, y=number, label="B", group=name), 
             color=SUCCESS_COLOR, fill='white')+
  xlab("Follow-up time (years)") + 
  ylab('Schedule') +
  scale_x_continuous(breaks=0:MAX_FOLLOW_UP,
                     limits = c(min_x, MAX_FOLLOW_UP)) +
  scale_y_continuous(limits = c(1 - 0.25, 4 + 0.25),
                     breaks = 1:4,
                     labels = schedule_df$name[!duplicated(schedule_df$name)],
                     sec.axis = dup_axis()) +
  theme(axis.title.x = element_blank(),
        axis.ticks.y.left = element_blank(),
        axis.title.y.left = element_blank(),
        axis.text.y.left = element_blank(),
        axis.text.y.right = element_text(size=FONT_SIZE),
        axis.title.y.right = element_text(size=FONT_SIZE))

##############
label_plot = ggplot() + 
  geom_vline(xintercept = c(0,last_biopsy_time,current_visit_time),
             color=c(WARNING_COLOR, SUCCESS_COLOR, "black"),
             linetype=c("solid", "solid", "dashed")) +
  geom_label(aes(x=c(0,last_biopsy_time,current_visit_time), 
                 y=c(0,0,0), 
                 label = c("Start\nsurveillance", 
                           "Last\nbiopsy",
                           "Current\nvisit")), color='white',
             size= LABEL_SIZE + 1.75,
             fill=c(WARNING_COLOR, SUCCESS_COLOR, 'black')) +
  theme_bw() +
  theme(text = element_text(size = FONT_SIZE),
        panel.background = element_blank(),
        axis.text = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        plot.margin = margin(0,0,0,0, unit = "pt")) + 
  xlab("Follow-up time (years)") + ylim(-0.1,0.1) + 
  xlim(min_x,MAX_FOLLOW_UP)

#######
expected_delays = sapply(split(schedule_df$times, schedule_df$name),
                         FUN = function(schedule){
                           set.seed(2019)
                           getConsequences(object = mvJoint_dre_psa_2knots_quad_age_light,
                                           patient_data = patient_df,
                                           cur_visit_time = cur_visit_time,
                                           latest_survival_time = last_biopsy_time,
                                           proposed_biopsy_times = schedule,
                                           M = 500, horizon = MAX_FOLLOW_UP)$expected_delay * 12
                         })

delay_plot = ggplot() + geom_col(aes(x=rep(1:4,2), 
                                     y=c(expected_delays, 12 - expected_delays)),
                                 color='black', fill=c(rep(c('darkgrey','white'), length(expected_delays))),
                                 width=0.5)+
  ylab("Expected time delay (months) in detecting progression") + 
  xlab("Schedule") +
  scale_y_continuous(breaks = seq(0, 12, length.out = 7),
                     labels = round(seq(0,12, length.out = 7),1),
                     limits=c(0,12))+
  scale_x_continuous(breaks=1:4, labels = names(expected_delays)) +
  coord_flip() +
  theme_bw() + 
  theme(text = element_text(size = FONT_SIZE))

demo_schedule = ggpubr::ggarrange(risk_plot, schedule_plot, 
                                  label_plot, delay_plot,
                                  ncol=1, align = "v",
                                  hjust = -1.5,
                                  heights = c(1,0.8,0.3,0.85),
                                  labels = c('A', 'B', '', 'C'), 
                                  legend = "none")

ggsave(demo_schedule, filename = "report/lastpaper/images/demo_schedule.eps",
       device = cairo_ps,  height = 7, width=7/1.1)
