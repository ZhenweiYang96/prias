mvJoint_dre_psa_dre_value$mcmc$b = NULL

DtList = seq(0.25, 10, 0.25)
thresholdsList = vector("list", length(DtList))
names(thresholdsList) = DtList

for(Dt in DtList){
  print(paste("Starting for Dt =", Dt))
  thresholdsList[[as.character(Dt)]] = find_thresholds.mvJMbayes_mod(mvJoint_dre_psa_dre_value, dre_psa_data_set, idVar = "P_ID", Dt = Dt, 
                                                                     n_cores = 4)
  print(paste("Done for Dt =", Dt))
  save(thresholdsList, file="Rdata/decision_analytic/DRE_PSA/thresholdsList.Rdata")
  print(paste("Saved for Dt =", Dt))
}
