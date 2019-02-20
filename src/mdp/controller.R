library(JMbayes)
library(survival)
library(splines)
library(ggplot2)

#Load a simulation
load("/home/a_tomer/Data/final_res_2nd_paper/jointModelData_seed_101_simNr_1_normal.Rdata")

fitted_JM = jointModelData$mvJoint_dre_psa_simDs
jointModelData$mvJoint_dre_psa_simDs = NULL

#choose a patient
patient_id = jointModelData$testData$testDs.id$P_ID[1]
patient_df = jointModelData$testData$testDs[jointModelData$testData$testDs$P_ID == patient_id,]
patient_df = patient_df[patient_df$visitTimeYears <=1,]

#source the common methods for all algorithms
source("src/mdp/common/simCommon.R")
#Source the method you want to use
source("src/mdp/tree_search/forward_search_no_Y.R")

#For DESPOT set these two
#N_DESPOT_SCENARIOS = 100
#DESPOT_TREE = list()

#We want to check how well we do with this approach. So we take a time point,
#and check how single threshold and forward search do for various patients


