library(JMbayes)
library(splines)
library(ggplot2)
library(ggpubr)

load("Rdata/gap3/PRIAS_2019/mvJoint_psa_time_scaled_light.Rdata")
load("Rdata/gap3/PRIAS_2019/cleandata.Rdata")
source("src/clinical_gap3/prediction_only_psa.R")

SUCCESS_COLOR = 'forestgreen'
DANGER_COLOR = 'red'
WARNING_COLOR = 'darkorange'
THEME_COLOR = 'dodgerblue4'
MAX_FOLLOW_UP = 10
POINT_SIZE = 2
FONT_SIZE = 12
LABEL_SIZE = 2.5

conditionalDynamicRiskPlot = function(object, pat_df, latest_survival_time=NA, 
                                      threshold = 0.1, xbreaks, xlabs, psa_breaks = NULL,
                                      max_follow_up){
  set.seed(2019)
  
  if(is.na(latest_survival_time)){
    latest_survival_time = max(pat_df$year_visit[!is.na(pat_df$gleason_sum)])
  }
  
  max_psa_time = max(pat_df$year_visit)
  print(max_psa_time)
  
  accuracy = 50
  psa_predict_times = seq(0, max_psa_time, length.out = accuracy)
  exp_fut = getExpectedFutureOutcomes(object, pat_df, latest_survival_time, Inf,
                                      psa_predict_times = psa_predict_times, 
                                      psaDist = "Tdist", addRandomError = F)
  mean_psa = rowMeans(exp_fut$predicted_psa)
  
  if(is.null(psa_breaks)){
    psa_breaks = seq(min(pat_df$log2psaplus1, na.rm = T), 
                     max(pat_df$log2psaplus1, na.rm = T), 
                     length.out = 3)
  }
  transformRiskToPSA = function(x){
    x*(tail(psa_breaks,1) - psa_breaks[1]) + psa_breaks[1]
  }
  
  #I am not expecting to show 10 tests but ok...its a reasonable upper limit for this code to not fail
  new_test_df = data.frame(survtime=numeric(), mean_cum_risk=numeric(),
                           mean_cum_risk_scaled=numeric(), 
                           lower_cum_risk_scaled=numeric(), upper_cum_risk_scaled=numeric())
  accuracy = 200
  latest_test_time = latest_survival_time
  while(latest_test_time<max_follow_up){
    survival_predict_times = seq(from = latest_test_time, to = max_follow_up, length.out = accuracy)
    exp_fut = getExpectedFutureOutcomes(object, pat_df, 
                                        latest_survival_time = latest_test_time, Inf,
                                        survival_predict_times = survival_predict_times, 
                                        psaDist = "Tdist", addRandomError = F)
    
    mean_cum_risk = 1-c(1, rowMeans(exp_fut$predicted_surv_prob))
    lower_cum_risk = 1-c(1, apply(exp_fut$predicted_surv_prob, 1, quantile, probs=0.025))
    upper_cum_risk = 1-c(1, apply(exp_fut$predicted_surv_prob, 1, quantile, probs=0.975))
    
    latest_test_time = survival_predict_times[which.min(abs(mean_cum_risk - threshold))]
    indices = survival_predict_times <= latest_test_time
    
    new_test_df = rbind(new_test_df, data.frame(survtime=survival_predict_times[indices],
                                                mean_cum_risk = mean_cum_risk[indices],
                                                mean_cum_risk_scaled = transformRiskToPSA(mean_cum_risk[indices]),
                                                lower_cum_risk_scaled = transformRiskToPSA(lower_cum_risk[indices]),
                                                upper_cum_risk_scaled = transformRiskToPSA(upper_cum_risk[indices])))
    
    print(latest_test_time)
  }
  
  riskBreaksOriginal = c(0, threshold, 0.5,1)
  riskAxisBreaks = transformRiskToPSA(riskBreaksOriginal)
  riskAxisLabels = paste0(riskBreaksOriginal*100, "%")
  riskAxisLabels[2] = paste0("Threshold\nk = ",riskAxisLabels[2])
  
  p = ggplot() +  
    geom_vline(xintercept = latest_survival_time, color=SUCCESS_COLOR)+
    geom_vline(xintercept = max_psa_time, linetype='dashed')+
    geom_segment(aes(x = latest_survival_time, xend=max_follow_up,
                     y=transformRiskToPSA(threshold), yend=transformRiskToPSA(threshold)),
                 linetype="dashed", color=DANGER_COLOR)+
    geom_line(aes(x=psa_predict_times, y=mean_psa), color=THEME_COLOR) + 
    geom_line(data=new_test_df, 
              aes(x=survtime, y=mean_cum_risk_scaled), color=DANGER_COLOR) +
    geom_ribbon(data=new_test_df, 
                aes(x=survtime, ymin=lower_cum_risk_scaled,
                    ymax=upper_cum_risk_scaled), alpha=0.15, fill=DANGER_COLOR) + 
    geom_vline(xintercept = new_test_df$survtime[new_test_df$mean_cum_risk==0],
               color=SUCCESS_COLOR) +
    geom_point(aes(x=pat_df$year_visit,y=pat_df$log2psaplus1),
               size=POINT_SIZE, color=THEME_COLOR) +
    theme_bw() + theme(text = element_text(size = FONT_SIZE), 
                       legend.title = element_blank(),
                       legend.position = "bottom",
                       axis.title.y = element_text(size=FONT_SIZE, color=THEME_COLOR),
                       axis.text.y = element_text(size=FONT_SIZE, color=THEME_COLOR),
                       axis.title.y.right  = element_text(size=FONT_SIZE, color=DANGER_COLOR),
                       axis.text.y.right = element_text(size=FONT_SIZE, color=DANGER_COLOR)) +
    scale_x_continuous(breaks = xbreaks, labels=xlabs,
                       limits = c(-0.35, max_follow_up), 
                       minor_breaks = seq(0, max_follow_up, 1)) +
    xlab("Follow-up time (years)")+
    scale_y_continuous(breaks = psa_breaks, 
                       labels = round(psa_breaks,1), 
                       limits = range(psa_breaks),
                       sec.axis = sec_axis(trans=~., 
                                           breaks= riskAxisBreaks,
                                           labels = riskAxisLabels,
                                           name = "Cumulative-risk of progression")) +
    ylab("Biomarker") +
    geom_point(aes(x=-5,y=-5, color="Observed PSA"), size=POINT_SIZE) +
    scale_color_manual(values = c(THEME_COLOR), labels="Observed longitudinal biomarker")
    
  
  return(p)
}

getpsaBreaks = function(pat_data, total=3){
  seq(min(pat_data$log2psaplus1, na.rm = T)-0.1,max(pat_data$log2psaplus1, na.rm = T) + 0.1, length.out = total)
}

set.seed(2019)
pat_data = prias_long_final[prias_long_final$P_ID==102,]
##I am perturbing the PSA of one of the patients to demo effect of rising PSA
pat_data$log2psaplus1[c(nrow(pat_data)-1, nrow(pat_data))] =  pat_data$log2psaplus1[c(nrow(pat_data)-1, nrow(pat_data))] + runif(n = 2, 0.1, 0.5)
psa_breaks = getpsaBreaks(pat_data)

cond_risk_plot = conditionalDynamicRiskPlot(mvJoint_psa_time_scaled,
                                            pat_data[pat_data$year_visit<=3.0,],
                                            latest_survival_time = 1.5,
                                            threshold = 0.15,
                                            xbreaks = 0:7,
                                            xlabs = 0:7,
                                            psa_breaks = psa_breaks,
                                            max_follow_up = 7) + 
  theme(axis.title.x = element_blank(),
        plot.margin = margin(0,0,0,0, unit = "pt"))

label_plot = ggplot() + 
  geom_vline(xintercept = c(0,1.5,2.589041,3.79397, 5.695033),
             color=c(WARNING_COLOR, SUCCESS_COLOR, "black",
                     SUCCESS_COLOR, SUCCESS_COLOR),
             linetype=c("solid", "solid", "dashed", "solid", "solid")) +
  geom_label(aes(x=c(0,1.5,2.589041), y=c(0,0,0), 
                 label = c("Start\nsurveillance", 
                           "Last\ntest",
                           "Current\nvisit")), color='white',
             size= LABEL_SIZE,
             fill=c(WARNING_COLOR, SUCCESS_COLOR, 'black')) +
  geom_label(aes(x=c(3.79397, 5.695033), y=c(0,0), 
                 label = c("1st future\ntest", 
                           "2nd future\ntest")), 
             color=SUCCESS_COLOR, size= LABEL_SIZE, fill='white') +
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
  xlim(-0.35,7)

schedule_plot = ggpubr::ggarrange(cond_risk_plot, label_plot,
                                  nrow=2, ncol=1, align="v",
                                  heights = c(4,1),
                                  common.legend = T, legend = "bottom")

print(schedule_plot)
ggsave(schedule_plot, filename = "report/lastpaper/images/schedule_explanation_102.eps",
       device = cairo_ps,  height = 4, width=6.5)
