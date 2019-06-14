library(JMbayes)
library(splines)
library(ggplot2)
library(ggpubr)
 
load("Rdata/gap3/PRIAS_2019/mvJoint_psa_time_scaled_light.Rdata")
load("Rdata/gap3/PRIAS_2019/cleandata.Rdata")
source("src/clinical_gap3/prediction_only_psa.R")

SUCCESS_COLOR = 'green'
DANGER_COLOR = 'red'
THEME_COLOR = 'dodgerblue4'
MAX_FOLLOW_UP = 10
POINT_SIZE = 2
FONT_SIZE = 15

jmExplanationPlot = function(object, pat_df, latest_survival_time=NA, 
                           xbreaks, xlabs,
                           max_follow_up){
  set.seed(2019)
  
  if(is.na(latest_survival_time)){
    latest_survival_time = max(pat_df$year_visit[!is.na(pat_df$gleason_sum)])
  }
  
  max_psa_time = max(pat_df$year_visit)
  
  accuracy = 50
  psa_predict_times = seq(0, max_psa_time, length.out = accuracy)
  survival_predict_times = seq(latest_survival_time, max_follow_up, length.out = accuracy)
  
  exp_fut = getExpectedFutureOutcomes(object, pat_df, latest_survival_time, Inf,
                                      survival_predict_times, psa_predict_times, psaDist = "Tdist")
  mean_psa = rowMeans(exp_fut$predicted_psa)
  mean_psa_velocity = rowMeans(exp_fut$predicted_psa_slope)
  mean_cum_risk = 1-c(1, rowMeans(exp_fut$predicted_surv_prob))
  lower_cum_risk = 1-c(1, apply(exp_fut$predicted_surv_prob, 1, quantile, probs=0.025))
  upper_cum_risk = 1-c(1, apply(exp_fut$predicted_surv_prob, 1, quantile, probs=0.975))
  
  common = ggplot() +  
    scale_x_continuous(breaks = xbreaks, labels=xlabs,
                       limits = c(0, max_follow_up), 
                       minor_breaks = seq(0, max_follow_up, 1)) +
    xlab("Follow-up time (years)")+
    geom_ribbon(aes(x=c(0,latest_survival_time), ymin=-c(Inf,Inf),
                    ymax=c(Inf, Inf)), alpha=0.15, fill=SUCCESS_COLOR) + 
    geom_vline(xintercept = latest_survival_time) + 
    geom_vline(xintercept = max_psa_time, linetype="dashed") +
    theme_bw() +
    theme(text = element_text(size = FONT_SIZE)) 
  
  blank_common = common + 
    theme(axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size=FONT_SIZE, color=THEME_COLOR),
          axis.text.y = element_text(size=FONT_SIZE, color=THEME_COLOR))
  
  A_y_breaks = seq(min(pat_df$log2psaplus1, na.rm = T), 
                   max(pat_df$log2psaplus1, na.rm = T), 
                   length.out = 3)
  
  A = blank_common + 
    geom_point(aes(x=pat_df$year_visit,y=pat_df$log2psaplus1),
               size=POINT_SIZE, color=THEME_COLOR) +
    geom_line(aes(x=psa_predict_times, y=mean_psa), color=THEME_COLOR) + 
    scale_y_continuous(breaks = A_y_breaks, labels = round(A_y_breaks,1)) +
    ylab(expression('log'[2]*'(PSA + 1) value'))
  
  B = blank_common + 
    geom_line(aes(x=psa_predict_times, y=mean_psa_velocity), color=THEME_COLOR) + 
    ylab(expression('log'[2]*'(PSA + 1) velocity')) 
  
  C = common + geom_line(aes(x=survival_predict_times, y=mean_cum_risk), color=DANGER_COLOR) +
    geom_ribbon(aes(x=survival_predict_times, ymin=lower_cum_risk,
                    ymax=upper_cum_risk), alpha=0.15, fill=DANGER_COLOR) + 
    geom_ribbon(aes(x=c(-1,-10),ymin=-Inf, ymax=Inf, fill="Region with Gleason ≤ 6"), alpha=0.15) +
    geom_point(aes(x=-5,y=-5, color="Observed PSA"), size=POINT_SIZE) +
    scale_fill_manual(values = c(SUCCESS_COLOR))+
    scale_color_manual(values = c(THEME_COLOR), labels=expression('Observed log'[2]*'(PSA + 1)'))+
    scale_y_continuous(breaks = seq(0,1, 0.5), 
                       labels = c("0%", "50%", "100%"), 
                       limits = c(0,1)) + ylab("Cumulative risk of\n Gleason \u2265 7") +
    theme(legend.title = element_blank(),
          axis.title.y = element_text(color=DANGER_COLOR, size=FONT_SIZE),
          axis.text.y = element_text(color=DANGER_COLOR, size=FONT_SIZE))
  
  upper_plot = ggarrange(A,B, C, ncol=1, nrow=3, align = "v",
                         heights = c(1,1,1.35), common.legend = T,
                        legend = "bottom", labels = "AUTO")
  return(upper_plot)
}

pat_data = prias_long_final[prias_long_final$P_ID==113,]
jm_explanation_plot = jmExplanationPlot(mvJoint_psa_time_scaled, pat_data,
                                latest_survival_time = 1, 
                                xbreaks = c(0, 1, 2, 3, pat_data$year_max_followup[1]),
                                xlabs = c("0\n(Start\nAS)", "1\n(Latest\nbiopsy)","2","3",
                                          "4\n (Latest\nvisit)"),
                                pat_data$year_max_followup[1])
ggsave(jm_explanation_plot, filename = "report/clinical/images/jmExplanationPlot_113.eps",
       device = cairo_ps, height = 9)