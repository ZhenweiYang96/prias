setwd("/home/a_tomer/Google Drive/PhD/src/prias/")
load("Rdata/decision_analytic/cleandata.Rdata")
source("src/decision_analytic/load_lib.R")

load("Rdata/decision_analytic/DRE_PSA/mvJoint_dre_psa_dre_value.Rdata")
mvJoint_dre_psa_dre_value$mcmc$b = NULL
dre_psa_data_set =  prias_long[!(is.na(prias_long$dre) & is.na(prias_long$psa)),]
dre_psa_data_set$high_dre = ifelse(dre_psa_data_set$dre=="T1c", 0, 1)
survfitJM(mvJoint_dre_psa_dre_value, dre_psa_data_set[dre_psa_data_set$P_ID %in% prias.id$P_ID[1],], idVar = "P_ID", survTimes = 10)
