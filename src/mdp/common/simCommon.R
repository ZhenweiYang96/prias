#Select parameters common for all types of algorithms
MIN_BIOPSY_GAP = 1

N_MCMC_ITER = 200
N_DESPOT_SCENARIOS = 100
DESPOT_TREE = list()

#Constants
MAX_FOLLOW_UP_TIME = 10
#Decision Epochs in years
PSA_CHECK_UP_TIME = c(seq(0, 2, 0.5), seq(2.5, MAX_FOLLOW_UP_TIME, 0.5), 15)
#DRE check up time years
DRE_CHECK_UP_TIME = c(seq(0, MAX_FOLLOW_UP_TIME, 0.5))

NR_DISCRETIZED_PSA = 3

LOWPSA_LOWDRE = "LP_LD"
MEDPSA_LOWDRE = "MP_LD"
HIGHPSA_LOWDRE = "HP_LD"
LOWPSA_HIGHDRE = "LP_HD"
MEDPSA_HIGHDRE = "MP_HD"
HIGHPSA_HIGHDRE = "HP_HD"

OUTCOME_CAT_NAMES = c(LOWPSA_LOWDRE, MEDPSA_LOWDRE, HIGHPSA_LOWDRE,
                      LOWPSA_HIGHDRE, MEDPSA_HIGHDRE, HIGHPSA_HIGHDRE)
OUTCOME_PSA_DRE_CAT = as.matrix(expand.grid(1:3, 0:1))
rownames(OUTCOME_PSA_DRE_CAT) = OUTCOME_CAT_NAMES

#Actions
BIOPSY = "B"
WAIT = "W"
ACTION_VECTOR = c(BIOPSY, WAIT)

#States
G6 = "G6"
G7 = "G7"
AT = "AT"
STATE_VECTOR = c(G6, G7, AT)

DISCOUNT_FACTOR = 1
DISCOUNT_FACTORS = DISCOUNT_FACTOR^((1:length(PSA_CHECK_UP_TIME)))
names(DISCOUNT_FACTORS) = PSA_CHECK_UP_TIME

TRUE_BIOPSY = "TB"
FALSE_BIOPSY = "FB"
TRUE_WAIT = "TW"
FALSE_WAIT = "FW"

reward_names = c(TRUE_BIOPSY, FALSE_BIOPSY, TRUE_WAIT, FALSE_WAIT)

REWARDS = c(5, -1, 1, -15)
names(REWARDS) = reward_names

source("src/mdp/common/prediction_psa_cat.R")

#########################################
# We define various functions from here onwards
#########################################
getNextDecisionEpoch = function(current_decision_epoch) {
  return(PSA_CHECK_UP_TIME[PSA_CHECK_UP_TIME>current_decision_epoch][1])
}

# getReward = function(current_state, action) {
#   if(current_state==AT){
#     return(0)
#   }else if(current_state==G7){
#     return(ifelse(action==BIOPSY, yes = REWARDS[TRUE_BIOPSY], no = REWARDS[FALSE_WAIT]))
#   }else{
#     return(ifelse(action==BIOPSY, yes = REWARDS[FALSE_BIOPSY], no = REWARDS[TRUE_WAIT]))
#   }
# }

# getReward = function(current_state, action, current_decision_epoch,
#                       latest_survival_time) {
#    if(current_state==AT){
#      return(0)
#    }else if(current_state==G7){
#      return(ifelse(action == BIOPSY, 
#                    yes = inv_time_to_biopsy_scale / (0.5 * (current_decision_epoch - latest_survival_time)),
#                    no = inv_time_to_biopsy_scale / (current_decision_epoch - getNextDecisionEpoch(current_decision_epoch))))
#    }else{
#      return(ifelse(action==BIOPSY, yes = -1, no = 1))
#    }
# }

getReward = function(current_state, action, current_decision_epoch,
                     latest_survival_time) {
  if(current_state==AT){
    return(0)
  }else if(current_state==G7){
    return(ifelse(action==BIOPSY,
                  yes = 0.5 * (current_decision_epoch - latest_survival_time) * time_to_biopsy_scale,
                  no = time_to_biopsy_scale * (current_decision_epoch - getNextDecisionEpoch(current_decision_epoch))))
  }else{
    return(ifelse(action==BIOPSY, yes = -1, no = 1))
  }
}