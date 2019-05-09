#Select parameters common for all types of algorithms
MIN_BIOPSY_GAP = 1

DESPOT_TREE = list()
DESPOT_BELIEF_CACHE = list()
DESPOT_Y_CACHE = list()

#Constants
MAX_FOLLOW_UP_TIME = 10
#Decision Epochs in years
PSA_CHECK_UP_TIME = c(seq(0, 2, 0.25), seq(2.5, MAX_FOLLOW_UP_TIME, 0.5), 15)
#DRE check up time years
DRE_CHECK_UP_TIME = c(seq(0, MAX_FOLLOW_UP_TIME, 0.5), 15)

BIOPSY_TEST_TIMES = DRE_CHECK_UP_TIME

NR_DISCRETIZED_PSA = 5

OUTCOME_CAT_NAMES = c(sapply(c("LD", "HD"), function(x){
  paste0(paste0("P", 1:NR_DISCRETIZED_PSA), x)
}))
OUTCOME_PSA_DRE_CAT = as.matrix(expand.grid(1:NR_DISCRETIZED_PSA, 0:1))
rownames(OUTCOME_PSA_DRE_CAT) = OUTCOME_CAT_NAMES

OUTCOME_TREATMENT = "OT"
OUTCOME_DUMMY = "OD"

#Actions
BIOPSY = "B"
WAIT = "W"
ACTION_VECTOR = c(BIOPSY, WAIT)

#States
G6 = "G6"
G7 = "G7"
AT = "AT"
STATE_VECTOR = c(G6, G7, AT)

TRUE_BIOPSY = "TB"
FALSE_BIOPSY = "FB"
TRUE_WAIT = "TW"
FALSE_WAIT = "FW"

reward_names = c(TRUE_BIOPSY, FALSE_BIOPSY, TRUE_WAIT, FALSE_WAIT)

#########################################
# We define various functions from here onwards
#########################################
getNextDecisionEpoch = function(current_decision_epoch) {
  return(BIOPSY_TEST_TIMES[BIOPSY_TEST_TIMES>current_decision_epoch][1])
}

getReward = function(current_state, action, ...) {
  if(current_state==AT){
    return(0)
  }else if(current_state==G7){
    return(ifelse(action==BIOPSY, yes = REWARDS[TRUE_BIOPSY], no = REWARDS[FALSE_WAIT]))
  }else{
    return(ifelse(action==BIOPSY, yes = REWARDS[FALSE_BIOPSY], no = REWARDS[TRUE_WAIT]))
  }
}

# getReward = function(current_state, action, current_decision_epoch,
#                      latest_survival_time, cur_biopsies) {
#   if(current_state==AT){
#     return(0)
#   }else if(current_state==G7){
#     return(ifelse(action==BIOPSY,
#                   yes = REWARDS[TRUE_BIOPSY],
#                   no = -(current_decision_epoch - latest_survival_time)))
#   }else{
#     return(ifelse(action==BIOPSY, 
#                   yes = cur_biopsies * REWARDS[FALSE_BIOPSY], 
#                   no = REWARDS[TRUE_WAIT]))
#   }
# }

getReward = function(current_state, action, current_decision_epoch,
                     latest_survival_time, cur_biopsies, conditionalFailTime) {
  if(current_state==AT){
    return(0)
  }else if(current_state==G7){
    
    delay = current_decision_epoch - conditionalFailTime
    
    #Adding negative reward of a biopsy to the true positive reward
    return(ifelse(action==BIOPSY,
                  yes = delay + REWARDS[FALSE_BIOPSY],
                  no = REWARDS[FALSE_WAIT] - delay))
    
  }else if(current_state==G6){
    
    #Assuming rewards of false biopsy is -1, and rewards of true wait is 0
    return(ifelse(action==BIOPSY,
                  yes = cur_biopsies * REWARDS[FALSE_BIOPSY],
                  no = REWARDS[TRUE_WAIT]))
    
  }
}

thresholdToReward = function(threshold, int_B, slope_B, int_W, slope_W){
  
  if(!is.na(threshold) & is.na(slope_W)){
    slope_W =  (int_B - int_W + slope_B * threshold) / (threshold)
  }
  
  R_FB = int_B
  R_TB = int_B + slope_B
  
  R_TW = int_W
  R_FW = int_W + slope_W
  
  rewards = c(R_TB, R_FB, R_TW, R_FW)
  names(rewards) = reward_names
  return(rewards)
}
