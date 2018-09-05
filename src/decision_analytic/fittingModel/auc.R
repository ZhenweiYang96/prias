auc1 = aucJM(mvJoint_dre_psa_dre_value, dre_psa_data_set, Tstart = 1, Thoriz=1 + 1.0101, idVar="P_ID")
save(auc1, file = "auc1.Rdata")
auc2 = aucJM(mvJoint_dre_psa_dre_value, dre_psa_data_set, Tstart = 2, Thoriz=2 + 1.0101, idVar="P_ID")
save(auc2, file = "auc2.Rdata")
auc3 = aucJM(mvJoint_dre_psa_dre_value, dre_psa_data_set, Tstart = 3, Thoriz=3 + 1.0101, idVar="P_ID")
save(auc3, file = "auc3.Rdata")


auc1 = aucJM(mvJoint_psa, psa_data_set, Tstart = 1, Thoriz=1 + 1.0101, idVar="P_ID")
save(auc1, file = "auc1_psa.Rdata")
auc2 = aucJM(mvJoint_psa, psa_data_set, Tstart = 2, Thoriz=2 + 1.0101, idVar="P_ID")
save(auc2, file = "auc2_psa.Rdata")
auc3 = aucJM(mvJoint_psa, psa_data_set, Tstart = 3, Thoriz=3 + 1.0101, idVar="P_ID")
save(auc3, file = "auc3_psa.Rdata")




