mvJoint_dre_psa_dre_value$mcmc$b = NULL

DtList = seq(0.5, 6, 0.5)
thresholdsList = vector("list", length(DtList))
names(thresholdsList) = DtList

for(Dt in DtList){
  thresholdsList[[as.character(Dt)]] = find_thresholds(mvJoint_dre_psa_dre_value, dre_psa_data_set, idVar = "P_ID", Dt = Dt, n_cores = 4)
  save(thresholdsList, file="Rdata/decision_analytic/DRE_PSA/thresholdsList.Rdata")
}
