library(JMbayes)
library(splines)
library(survival)
library(interval)

load("Rdata/gap3/PRIAS_2019/validation/predicted_risk_comparisons/UCSF.Rdata")
load("Rdata/gap3/PRIAS_2019/npmle_all.Rdata")

npmle=npmle_all$`London-KCL`
npmle_time_points = as.numeric(cbind(c(0,0), npmle$intmap))
npmle_cumrisk = c(0, as.numeric(rep(c(0, cumsum(npmle$pf)), each=2)))[1:length(npmle_time_points)]

cohort_model_pred = cumrisk_models[[1]]
prias_model_pred = cumrisk_models[[2]]
prias_model_recalib_pred = cumrisk_models[[3]]

calib_pred_times = seq(0, 10, 0.1)

rm(list = setdiff(ls(), c("cohort_model_pred", "prias_model_pred", "prias_model_recalib_pred",
                          "calib_pred_times", "npmle_time_points", "npmle_cumrisk")))

#Checking the Kaplan Meier
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
  ylab("Cumulative-risk of Reclassification") + 
  xlab("Follow-up time (years)") +
  theme_bw() + 
  theme(legend.position = "bottom", 
        legend.title = element_blank(),
        text=element_text(size=15))
A = .Last.value
#Checking the individual predictions
ggplot() + 
  geom_boxplot(aes(x=rep("PRIAS", prod(dim(prias_model_pred))), 
                   y=c(prias_model_pred-cohort_model_pred)), outlier.shape = NA) +
  geom_boxplot(aes(x=rep("PRIAS Recalib", prod(dim(prias_model_recalib_pred))), 
                   y=c(prias_model_recalib_pred-cohort_model_pred)), outlier.shape = NA) +
  scale_y_continuous(breaks = seq(-1,1,0.2), limits = c(-1,1)) +
  geom_hline(yintercept = 0, color='red') +
  ylab("Prediction Error") + 
  xlab("Model") +
  theme_bw() + 
  theme(legend.position = "bottom", 
        legend.title = element_blank(),
        text=element_text(size=15))
B = .Last.value

ggpubr::ggarrange(A,B, ncol=2, nrow=1)
