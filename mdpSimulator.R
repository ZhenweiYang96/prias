#Libraries
library(splines)
library(survival)
library(JMbayes)

#Rdata
load("~/Data/final_res_2nd_paper/jointModelData_seed_101_simNr_1_normal.Rdata")

#source code
source("src/mdp/postRandEff.R")
source("src/mdp/mcts.R")
source("src/mdp/survfitJM_noSeed.R")

fitted_JM = jointModelData$mvJoint_dre_psa_simDs
patient_data = jointModelData$testData$testDs[jointModelData$testData$testDs$P_ID==893,]
patient_data_starting = patient_data[patient_data$visitTimeYears<=1,]

#Set some constants
#This one is currently not used
EXPLORATION_REWARD_CONSTANT = 100

DISCOUNT_FACTOR = 1

tt = selectAction(patient_obs_df = patient_data_starting, 
             n_simulations = 50000, current_decision_epoch = 1, latest_survival_time = 0,
             max_decision_epoch = 2.5)

tt_noobs = selectAction(patient_obs_df = patient_data_starting, 
                  n_simulations = 50000, current_decision_epoch = 1, latest_survival_time = 0,
                  max_decision_epoch = 2.5)
