model_sqrt = mvglmer(list(
  sqrtpsa ~  I(Age - 70) +  I((Age - 70)^2) + 
    ns(visitTimeYears, knots=c(0.1, 0.5, 4), Boundary.knots=c(0, 7)) + 
    (ns(visitTimeYears, knots=c(0.1), Boundary.knots=c(0, 7))|P_ID)), 
  data=training_psa_data_set, families = list(gaussian))

save(model_sqrt, file="Rdata/Gleason as event/tdist/modelsqrt.Rdata")


model = mvglm_psa_spline_pt1pt54_pt1

fittedmarginal = model$components$X1 %*% model$postMeans$betas1
fittedSubjects = rowSums(model$components$Z1 * model$postMeans$b[as.numeric(droplevels(training_psa_data_set$P_ID)),,drop=FALSE])

training_psa_data_set$fitted = fittedmarginal + fittedSubjects
training_psa_data_set$resid = training_psa_data_set$log2psa - training_psa_data_set$fitted

training_psa_data_set$standardized_resid = (training_psa_data_set$resid - mean(training_psa_data_set$resid))/model$postMeans$sigma1  

#Both qqplots give same result

#QQPLOT after explicitly standardizing the residuals
ggplot(data=training_psa_data_set) + geom_qq(aes(sample=standardized_resid), dparams = list(df=3),distribution = qt) + geom_abline(intercept = 0, slope = 1)

#QQPLOT after scaling theoretical distribution
ggplot(data=training_psa_data_set) + geom_qq(aes(sample=resid), 
                                             dparams = list(mu=0, sigma=model$postMeans$sigma1, df=3),
                                             distribution = JMbayes:::qgt) + 
  xlim(c(-5,5)) + 
  ylim(c(-5,5)) +
  geom_abline(intercept = -0.0636430822997942, slope = 0.843009408103346) + 
  xlab("Theoretical quantiles") + ylab("Residual quantiles")

p2 = ggplot(data=training_psa_data_set) + geom_qq(aes(sample=resid), 
                                             dparams = list(df=3),
                                             distribution = qt) + 
   geom_abline(intercept = 0.00153022619620519, slope = 0.152554094219411) + 
  xlab("T-distribution (df=3) quantiles") + ylab("Residual quantiles") + 
  ggtitle("T-distributed (df=3) errors") + 
  theme(text = element_text(size=11), axis.text=element_text(size=11), 
        plot.title = element_text(hjust = 0.5, size=13)) 

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