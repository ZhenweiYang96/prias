#Select parameters common for all types of algorithms
MIN_BIOPSY_GAP = 1

#Constants
MAX_FOLLOW_UP_TIME = 13.5
#Decision Epochs in years
PSA_CHECK_UP_TIME = c(seq(0, 2, 0.25), seq(2.5, MAX_FOLLOW_UP_TIME, 0.5))
#DRE check up time years
DRE_CHECK_UP_TIME = seq(0, MAX_FOLLOW_UP_TIME, 0.5)

BIOPSY_TEST_TIMES = DRE_CHECK_UP_TIME

getNextDecisionEpoch = function(current_decision_epoch) {
  return(BIOPSY_TEST_TIMES[BIOPSY_TEST_TIMES>current_decision_epoch][1])
}
