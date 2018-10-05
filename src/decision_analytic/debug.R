#Step 1: AUC and Prederr on Real prias data with real prias fitted model object
auc_0_0pt5 =  aucJM.mvJMbayes_mod(mvJoint_dre_psa_dre_value, newdata = prias_long, 
                              Tstart = 0.0001, Thoriz = 0.5001, idVar = "P_ID")
auc_1_1pt5 = aucJM.mvJMbayes_mod(mvJoint_dre_psa_dre_value, newdata = prias_long, 
                              Tstart = 1.0001, Thoriz = 1.5001, idVar = "P_ID")
auc_2_2pt5 = aucJM.mvJMbayes_mod(mvJoint_dre_psa_dre_value, newdata = prias_long, 
                              Tstart = 2.0001, Thoriz = 2.5001, idVar = "P_ID")
auc_3_3pt5 = aucJM.mvJMbayes_mod(mvJoint_dre_psa_dre_value, newdata = prias_long, 
                              Tstart = 3.0001, Thoriz = 3.5001, idVar = "P_ID")
auc_4_4pt5 = aucJM.mvJMbayes_mod(mvJoint_dre_psa_dre_value, newdata = prias_long, 
                              Tstart = 4.0001, Thoriz = 4.5001, idVar = "P_ID")

save.image(file="Rdata/decision_analytic/AUC_Brier/auc_prederr_t_pt5.Rdata")

prederr_0_0pt5 =  prederrJM.mvJMbayes_mod(mvJoint_dre_psa_dre_value, newdata = prias_long, 
                                      Tstart = 0.0001, Thoriz = 0.5001, idVar = "P_ID")
prederr_1_1pt5 = prederrJM.mvJMbayes_mod(mvJoint_dre_psa_dre_value, newdata = prias_long, 
                                      Tstart = 1.0001, Thoriz = 1.5001, idVar = "P_ID")
prederr_2_2pt5 = prederrJM.mvJMbayes_mod(mvJoint_dre_psa_dre_value, newdata = prias_long, 
                                      Tstart = 2.0001, Thoriz = 2.5001, idVar = "P_ID")
prederr_3_3pt5 = prederrJM.mvJMbayes_mod(mvJoint_dre_psa_dre_value, newdata = prias_long, 
                                      Tstart = 3.0001, Thoriz = 3.5001, idVar = "P_ID")
prederr_4_4pt5 = prederrJM.mvJMbayes_mod(mvJoint_dre_psa_dre_value, newdata = prias_long, 
                                      Tstart = 4.0001, Thoriz = 4.5001, idVar = "P_ID")

save.image(file="Rdata/decision_analytic/AUC_Brier/auc_prederr_t_pt5.Rdata")

