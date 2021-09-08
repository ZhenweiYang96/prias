library(JMbayes)
library(splines)
library(ggplot2)
library(ggpubr)

load("Rdata/gap3/PRIAS_2019/mvJoint_psa_time_scaled_light.Rdata")
load("Rdata/gap3/PRIAS_2019/cleandata.Rdata")
source("src/clinical_gap3/prediction_only_psa.R")

SUCCESS_COLOR = 'forestgreen'
DANGER_COLOR = 'red'
THEME_COLOR = 'dodgerblue4'
MAX_FOLLOW_UP = 10
POINT_SIZE = 2
FONT_SIZE = 12

dynamicRiskPlot = function(object, pat_df, latest_survival_time=NA, 
                           xbreaks, xlabs, psa_breaks = NULL,
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
  
  if(is.null(psa_breaks)){
    psa_breaks = seq(min(pat_df$log2psaplus1, na.rm = T), 
                     max(pat_df$log2psaplus1, na.rm = T), 
                     length.out = 3)
  }
  transformRiskToPSA = function(x){
    x*(tail(psa_breaks,1) - psa_breaks[1]) + psa_breaks[1]
  }
  
  mean_cum_risk = transformRiskToPSA(1-c(1, rowMeans(exp_fut$predicted_surv_prob)))
  lower_cum_risk = transformRiskToPSA(1-c(1, apply(exp_fut$predicted_surv_prob, 1, quantile, probs=0.025)))
  upper_cum_risk = transformRiskToPSA(1-c(1, apply(exp_fut$predicted_surv_prob, 1, quantile, probs=0.975)))
  
  riskAxisBreaks = transformRiskToPSA(c(0, 0.5, 1))
  riskAxisLabels = c("0%", "50%", "100%")
  
  p = ggplot() +  
    scale_x_continuous(breaks = xbreaks, labels=xlabs,
                       limits = c(0, max_follow_up), 
                       minor_breaks = seq(0, max_follow_up, 1)) +
    xlab("Follow-up time (years)")+
    geom_vline(xintercept = latest_survival_time, color=SUCCESS_COLOR) + 
    geom_vline(xintercept = max_psa_time, linetype="dashed") +
    geom_point(aes(x=pat_df$year_visit,y=pat_df$log2psaplus1),
               size=POINT_SIZE, color=THEME_COLOR) +
    geom_line(aes(x=psa_predict_times, y=mean_psa), color=THEME_COLOR) + 
    geom_line(aes(x=survival_predict_times, y=mean_cum_risk), color=DANGER_COLOR, linetype='dotted') +
    geom_ribbon(aes(x=survival_predict_times, ymin=lower_cum_risk,
                    ymax=upper_cum_risk), alpha=0.15, fill=DANGER_COLOR) + 
    scale_y_continuous(breaks = psa_breaks, 
                       labels = round(psa_breaks,1), 
                       limits = range(psa_breaks),
                       sec.axis = sec_axis(trans=~., 
                                           breaks= riskAxisBreaks,
                                           labels = riskAxisLabels,
                                           name = "Cumulative-risk of\n progression")) +
    geom_point(aes(x=-5,y=-5, color="Observed PSA"), size=POINT_SIZE) +
    scale_color_manual(values = c(THEME_COLOR), labels="Observed longitudinal biomarker")+
    ylab("Biomarker") + 
    theme_bw() + theme(text = element_text(size = FONT_SIZE), 
                       legend.title = element_blank(),
                       axis.title.y = element_text(size=FONT_SIZE, color=THEME_COLOR),
                       axis.text.y = element_text(size=FONT_SIZE, color=THEME_COLOR),
                       axis.title.y.right  = element_text(size=FONT_SIZE, color=DANGER_COLOR),
                       axis.text.y.right = element_text(size=FONT_SIZE, color=DANGER_COLOR))
  
  return(p)
}

getpsaBreaks = function(pat_data, total=3){
  seq(min(pat_data$log2psaplus1, na.rm = T)-0.1,max(pat_data$log2psaplus1, na.rm = T) + 0.1, length.out = total)
}

set.seed(2019)
pat_data = prias_long_final[prias_long_final$P_ID==102,]
pat_data = pat_data[1:16,]
##I am perturbing the PSA of one of the patients to demo effect of rising PSA
pat_data$log2psaplus1[c(nrow(pat_data)-1, nrow(pat_data))] =  pat_data$log2psaplus1[c(nrow(pat_data)-1, nrow(pat_data))] + runif(n = 2, 0.1, 0.5)
pat_data$year_visit[8]=2.5
pat_data$year_visit[10]=3.5
pat_data$year_visit[16]=6.5
psa_breaks = getpsaBreaks(pat_data)

dynrisk_plot_1 = dynamicRiskPlot(mvJoint_psa_time_scaled,
                                 pat_data[pat_data$year_visit<=2.5,],
                                 latest_survival_time = 1.5,
                                 xbreaks = c(0, 1.5, 2.5, 6.5),
                                 xlabs = c("0\n(Start\nSurveillance)", "t=1.5\n(Last\ntest)",
                                           "v=2.5\n(Current\nvisit)","6.5"),
                                 psa_breaks = psa_breaks,
                                 max_follow_up = 6.5) + theme(axis.title.x = element_blank())

dynrisk_plot_2 = dynamicRiskPlot(mvJoint_psa_time_scaled,
                                 pat_data[pat_data$year_visit<=3.5,],
                                 latest_survival_time = 2.5,
                                 xbreaks = c(0, 2.5, 3.5, 6.5),
                                 xlabs = c("0\n(Start\nSurveillance)", "t=2.5\n (Last\ntest)",
                                           "v=3.5\n(Current\nvisit)","6.5"),
                                 psa_breaks = psa_breaks,
                                 max_follow_up = 6.5) + theme(axis.title.x = element_blank())

dynrisk_plot_3 = dynamicRiskPlot(mvJoint_psa_time_scaled,
                                 pat_data,
                                 latest_survival_time = 2.5,
                                 xbreaks = c(0, 2.5, 6.5),
                                 xlabs = c("0\n(Start\nSurveillance)", "t=2.5\n (Last\ntest)",
                                           "v=6.5\n (Current\nvisit)"),
                                 psa_breaks = psa_breaks,
                                 max_follow_up = 6.5)

dynrisk_plot = ggarrange(dynrisk_plot_1, dynrisk_plot_2, dynrisk_plot_3,
          align = "v", labels = "AUTO", heights = c(1,1,1.1),
          nrow = 3, ncol=1, legend = "bottom", common.legend = T)

print(dynrisk_plot) 
ggsave(dynrisk_plot, filename = "report/lastpaper/images/dynrisk_plot_102.pdf",
       device = cairo_pdf,  height = 7, width=7/1.333)
ggsave(dynrisk_plot, filename = "report/lastpaper/figure2.eps",
       device = cairo_ps,  height = 7, width=7/1.333)