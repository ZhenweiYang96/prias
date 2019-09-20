library(JMbayes)
library(splines)
library(survival)
library(interval)

load("Rdata/gap3/PRIAS_2019/validation/predicted_risk_comparisons/Toronto.Rdata")
load("Rdata/gap3/PRIAS_2019/npmle_all.Rdata")

npmle=npmle_all$Toronto
npmle_time_points = as.numeric(cbind(c(0,0), npmle$intmap))
npmle_cumrisk = c(0, as.numeric(rep(c(0, cumsum(npmle$pf)), each=2)))[1:length(npmle_time_points)]

cohort_model_pred = cumrisk_models[[1]]
prias_model_pred = cumrisk_models[[2]]
prias_model_recalib_pred = cumrisk_models[[3]]

calib_pred_times = seq(0, 10, 0.1)

rm(list = setdiff(ls(), c("cohort_model_pred", "prias_model_pred", "prias_model_recalib_pred",
                          "calib_pred_times", "npmle_time_points", "npmle_cumrisk")))

ggplot() + 
  geom_line(aes(x=npmle_time_points, y=npmle_cumrisk, color="NPMLE")) + 
  geom_line(aes(x=calib_pred_times, 
                y=rowMeans(cohort_model_pred, na.rm = T), color="Cohort Model")) +
  geom_line(aes(x=calib_pred_times, 
                y=rowMeans(prias_model_pred, na.rm = T), color="PRIAS Model")) + 
  geom_line(aes(x=calib_pred_times, 
                y=rowMeans(prias_model_recalib_pred, na.rm = T), color="PRIAS Recalib Model")) + 
  scale_x_continuous(breaks=0:10, limits=c(0,10)) + 
  scale_y_continuous(breaks = seq(0,1, 0.25), labels=paste0(seq(0,1,0.25)*100, "%")) + 
  ylab("Cumulative-risk of Gleason > 6") + 
  xlab("Follow-up time (years)") +
  theme_bw() + 
  theme(legend.position = "bottom", 
        legend.title = element_blank(),
        text=element_text(size=15))
