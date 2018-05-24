model_sqrt = mvglmer(list(
  sqrtpsa ~  I(Age - 70) +  I((Age - 70)^2) + 
    ns(visitTimeYears, knots=c(0.1, 0.5, 4), Boundary.knots=c(0, 7)) + 
    (ns(visitTimeYears, knots=c(0.1), Boundary.knots=c(0, 7))|P_ID)), 
  data=training_psa_data_set, families = list(gaussian))

save(model_sqrt, file="Rdata/Gleason as event/tdist/modelsqrt.Rdata")

modelList = list(mvglm_psa_spline_pt1pt54_pt1_t3, mvglm_log2psa_pluspt1_spline_pt1pt54_pt1,
                 mvglm_log2psa_plus1_spline_pt1pt54_pt1)

resid_psa = lapply(modelList, function(model){
  fittedmarginal = model$components$X1 %*% model$postMeans$betas1
  fittedSubjects = rowSums(model$components$Z1 * model$postMeans$b[as.numeric(droplevels(model$data$P_ID)),,drop=FALSE])
  fitted_psa = fittedmarginal + fittedSubjects
  resid_psa = model$components$y1 - fitted_psa
  resid_psa
})

source("src/pers_sched/qqtdist.R")

lapply(resid_psa, function(residual_psa){
  qqtdist(residual_psa, df=3)
  qqlineTdist(residual_psa, df=3)
})

p1 = ggplot() + geom_qq(aes(sample=resid_psa[[1]]), 
                                             dparams = list(df=3),
                                             distribution = qt) + 
   geom_abline(intercept = 0.00153022619620519, slope = 0.152554094219411) + 
  xlab("T-distribution (df=3) quantiles") + ylab("Residual quantiles") + 
  ggtitle(expression('Using '*'log'[2]*' (PSA)'*' transformation')) + 
  theme(text = element_text(size=11), axis.text=element_text(size=11), 
        plot.title = element_text(hjust = 0.5, size=13)) 
print(p1)

p2 = ggplot() + geom_qq(aes(sample=resid_psa[[2]]), 
                        dparams = list(df=3),
                        distribution = qt) + 
  geom_abline(intercept = 0.00131744150956975, slope = 0.149009696633965) + 
  xlab("T-distribution (df=3) quantiles") + ylab("Residual quantiles") + 
  ggtitle(expression('Using '*'log'[2]*' (PSA + 0.1)'*' transformation')) + 
  theme(text = element_text(size=11), axis.text=element_text(size=11), 
        plot.title = element_text(hjust = 0.5, size=13)) 
print(p2)

p3 = ggplot() + geom_qq(aes(sample=resid_psa[[3]]), 
                        dparams = list(df=3),
                        distribution = qt) + 
  geom_abline(intercept = 0.00159500079649591, slope = 0.125110597597315) + 
  xlab("T-distribution (df=3) quantiles") + ylab("Residual quantiles") + 
  ggtitle(expression('Using '*'log'[2]*' (PSA + 1)'*' transformation')) + 
  theme(text = element_text(size=11), axis.text=element_text(size=11), 
        plot.title = element_text(hjust = 0.5, size=13)) 
print(p3)

ggsave(file="report/pers_sched/latex/biometrics_submission/images/model_fit/qqplot_various_log_transform_t3.eps", 
       multiplot(p1,p3, p2, cols=2),
       width=8.27, height=8.27)

source("../JMBayes/Anirudh/dev/multiplot.R")

ggsave(file="report/pers_schedule/biometrics_submission/images/qqplot_t3.eps", 
       width=8.27, height=9.69/1.25, device=cairo_ps)

ggsave(file="report/pers_schedule/biometrics_submission/images/qqplot_norm_t3.eps", 
       multiplot(p1,p2, cols=2),
       width=8.27, height=8.27/2)

###################
joint_psa_replaced_prias_t3$Funs$densLong = function(y, eta.y, scale, log = FALSE, data){dgt(x=y, mu=eta.y, sigma=scale, df = 3, log = log)}

##################
alphaDf = data.frame(type=c("Value", "Slope"), mean=mvJoint_psa_spline_pt1pt54_pt1_tdboth_t3$statistics$postMeans$alphas, low=mvJoint_psa_spline_pt1pt54_pt1_tdboth_t3$statistics$CIs$alphas[1,], up=mvJoint_psa_spline_pt1pt54_pt1_tdboth_t3$statistics$CIs$alphas[2,])
ggplot(data=alphaDf) + geom_point(aes(type, mean), size=4) + 
  geom_errorbar(width=.1, aes(type, ymin=low, ymax=up), size=0.9) + 
  ylab("log Hazard Ratio") + xlab("Association Parameter") + 
  geom_hline(yintercept = 0, linetype="dashed") +
  theme(text = element_text(size=25), axis.text=element_text(size=25))  