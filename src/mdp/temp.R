environment(get_b_fullBayes) = asNamespace("JMbayes")
environment(predictLongitudinalOutcome) = asNamespace("JMbayes")

mse_1 = c()
mse_2 = c()
for(P_ID in unique(dre_psa_data_set$P_ID)[1:100]){
  #print(P_ID)
  patient_132 = dre_psa_data_set[dre_psa_data_set$P_ID==P_ID,]
  p_time_start = patient_132$progression_time_start[1]
  p_time_end = patient_132$progression_time_end[1]
  
  greg_b = survfit_temp(mvJoint_dre_psa_dre_value, patient_132, last.time = p_time_end, idVar = "P_ID")
  t_start = Sys.time()
  optim_res = get_b_fullBayes(object = mvJoint_dre_psa_dre_value,
                              patient_data = patient_132,
                              timeLowerLimit = p_time_end,
                              timeUpperLimit = Inf, M = 0)
  b_1 = optim_res$par
  
  t_end = Sys.time()
  t_end-t_start
  #print(b_1)
  # b_1 = getPostRandEff(object = mvJoint_dre_psa_dre_value, 
  #                      patient_data = patient_132,
  #                      timeLowerLimit = 0,
  #                      timeUpperLimit = Inf)
  
  ggplot() + geom_point(aes(x=patient_132$visitTimeYears, y=patient_132$log2psaplus1)) +
     geom_line(aes(x=patient_132$visitTimeYears, y=psaXbetaZb(patient_132$Age[1],patient_132$visitTimeYears,mvJoint_dre_psa_dre_value$statistics$postMeans$betas2, b_1[3:7]), color="b1"))+
    geom_line(aes(x=patient_132$visitTimeYears, y=psaXbetaZb(patient_132$Age[1],patient_132$visitTimeYears,mvJoint_dre_psa_dre_value$statistics$postMeans$betas2, greg_b[3:7]), color="b2"))
  
  mse_1=c(mse_1, patient_132$log2psaplus1 - psaXbetaZb(patient_132$Age[1],patient_132$visitTimeYears,mvJoint_dre_psa_dre_value$statistics$postMeans$betas2, b_1[3:7]))
  mse_2=c(mse_2, patient_132$log2psaplus1 - psaXbetaZb(patient_132$Age[1],patient_132$visitTimeYears,mvJoint_dre_psa_dre_value$statistics$postMeans$betas2, greg_b[3:7]))
}
ggplot() + geom_point(aes(x=patient_132$visitTimeYears, y=patient_132$high_dre)) +
  geom_line(aes(x=patient_132$visitTimeYears, y=plogis(dreLogOddsXbetaZb(patient_132$Age[1],patient_132$visitTimeYears,mvJoint_dre_psa_dre_value$statistics$postMeans$betas1, b_1[1:2])), color="b1"))+
  geom_line(aes(x=patient_132$visitTimeYears, y=plogis(dreLogOddsXbetaZb(patient_132$Age[1],patient_132$visitTimeYears,mvJoint_dre_psa_dre_value$statistics$postMeans$betas1, greg_b[1:2])), color="b2"))

t_start = Sys.time()
fitRes = predictLongitudinalOutcome(object = mvJoint_dre_psa_dre_value,
                                    patient_data = patient_132,
                                    timeLowerLimit = p_time_start,
                                    timeUpperLimit = p_time_end, M = 500, visitTimeYears = patient_132$visitTimeYears)
t_end = Sys.time()
t_end-t_start

ggplot() + geom_point(aes(x=patient_132$visitTimeYears, y=patient_132$log2psaplus1)) + 
  geom_line(aes(x=patient_132$visitTimeYears, y=apply(fitRes$predicted_psa, 1, mean))) +
  geom_ribbon(aes(x=patient_132$visitTimeYears, 
                  ymin=apply(fitRes$predicted_psa, 1, quantile, probs=0.025), 
                  ymax=apply(fitRes$predicted_psa, 1, quantile, probs=0.975)),
              fill="gray", alpha=0.5) + xlim(0,7) + ylim(1.8,4)

ggplot() + geom_point(aes(x=patient_132$visitTimeYears, y=patient_132$high_dre)) + 
  geom_line(aes(x=patient_132$visitTimeYears, y=apply(plogis(fitRes$predicted_dreLogOdds), 1, mean))) +
  geom_ribbon(aes(x=patient_132$visitTimeYears, 
                  ymin=apply(plogis(fitRes$predicted_dreLogOdds), 1, quantile, probs=0.025), 
                  ymax=apply(plogis(fitRes$predicted_dreLogOdds), 1, quantile, probs=0.975)),
              fill="gray", alpha=0.5) + xlim(0,7) 

ggplot() + geom_point(aes(x=patient_132$visitTimeYears, y=patient_132$high_dre)) + 
  geom_line(aes(x=patient_132$visitTimeYears, y=apply(plogis(fitRes$predicted_dreLogOdds), 1, mean), color="A")) +
  geom_line(aes(x=greg_fitRes$visitTimeYears, y=apply(plogis(greg_fitRes$trueDRELogOdds), 1, mean), color="G")) +
  xlim(0,7) 

ggplot() + geom_point(aes(x=patient_132$visitTimeYears, y=patient_132$log2psaplus1)) + 
  geom_line(aes(x=patient_132$visitTimeYears, y=apply((fitRes$predicted_psa), 1, mean), color="A")) +
  geom_line(aes(x=greg_fitRes$visitTimeYears, y=apply((greg_fitRes$trueLog2psaplus1), 1, mean), color="G")) +
  xlim(0,7) + ylim(1.8, 4)

t_start = Sys.time()
greg_fitRes = predictPSADRE(object = mvJoint_dre_psa_dre_value,
                                    newdata = patient_132,
                                    last.time = p_time_end, idVar = "P_ID",
                                    M = 500, seed = 2019)
t_end = Sys.time()
t_end-t_start


ggplot() + geom_point(aes(x=patient_132$visitTimeYears, y=patient_132$log2psaplus1)) + 
  geom_line(aes(x=greg_fitRes$visitTimeYears, y=apply(greg_fitRes$trueLog2psaplus1, 1, mean))) +
  geom_ribbon(aes(x=greg_fitRes$visitTimeYears, 
                  ymin=apply(greg_fitRes$trueLog2psaplus1, 1, quantile, probs=0.025), 
                  ymax=apply(greg_fitRes$trueLog2psaplus1, 1, quantile, probs=0.975)),
              fill="gray", alpha=0.5) + xlim(0,7) + ylim(1.8,4)

ggplot() + geom_point(aes(x=patient_132$visitTimeYears, y=patient_132$high_dre)) + 
  geom_line(aes(x=greg_fitRes$visitTimeYears, y=apply(plogis(greg_fitRes$trueDRELogOdds), 1, mean))) +
  geom_ribbon(aes(x=greg_fitRes$visitTimeYears, 
                  ymin=apply(plogis(greg_fitRes$trueDRELogOdds), 1, quantile, probs=0.025), 
                  ymax=apply(plogis(greg_fitRes$trueDRELogOdds), 1, quantile, probs=0.975)),
              fill="gray", alpha=0.5) + xlim(0,7) 
