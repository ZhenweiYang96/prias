library(JMbayes)
library(survival)
library(splines)
library(ggplot2)

#source the common methods for all algorithms
source("src/mdp/common/simCommon.R")
source("src/mdp/common/prediction.R")

max_cores = 3
N_MCMC_ITER = 500

file = list.files(path='C://Users//838035//prias//Rdata//mdp//final_res', full.names = T)[2]
load(file)
print(paste("Loaded file:", 1))

fitted_JM = jointModelData$mvJoint_dre_psa_simDs
testData = jointModelData$testData
jointModelData$mvJoint_dre_psa_simDs = NULL

threshold = 0.1

ct = makeCluster(max_cores)
registerDoParallel(ct)

res = foreach(pid=testData$testDs.id$P_ID,
              .packages = c("splines", "JMbayes"), .combine = 'rbind') %dopar%{
                source("src/mdp/common/prediction.R")
                
                set.seed(pid)
                
                patient_df = testData$testDs[testData$testDs$P_ID==pid,]
                progression_time = patient_df$progression_time[1]
                
                decision_epoch = 2
                pat_data = patient_df[patient_df$visitTimeYears <= decision_epoch,]
                
                survTimes = 1:10
                survProbMatrix = getExpectedFutureOutcomes(fitted_JM, pat_data, 0, Inf, survTimes, M = N_MCMC_ITER)$predicted_surv_prob
                
                meanSurvProb = rowMeans(survProbMatrix)
                lowSurvProb = apply(survProbMatrix, 1, quantile, probs=0.025)
                uppSurvProb = apply(survProbMatrix, 1, quantile, probs=0.975)
                
                wtPoints = getGaussianQuadWeightsPoints(c(0,10))
                points = sort(wtPoints$points, decreasing = F)
                wt = wtPoints$weights[order(wtPoints$points, decreasing = F)]
                
                survProbMatrix = getExpectedFutureOutcomes(fitted_JM, pat_data, 0, Inf, points, M = N_MCMC_ITER)$predicted_surv_prob
                auc = sum(wt*rowMeans(survProbMatrix))
              
                condSurvTimes = 3:10
                condSurvProbMatrix = getExpectedFutureOutcomes(fitted_JM, pat_data, latest_survival_time = 2, earliest_failure_time = Inf, 
                                                               condSurvTimes, M = N_MCMC_ITER)$predicted_surv_prob
                
                meanCondSurvProb = c(1,1, rowMeans(condSurvProbMatrix))
                
                return(data.frame('P_ID'=pid,'progression_time'=progression_time,
                                  'auc'=auc, 'survTimes'=survTimes, 'meanSurvProb'=meanSurvProb,
                                  'lowSurvProb'=lowSurvProb, 'uppSurvProb'=uppSurvProb,
                                  'meanCondSurvProb'=meanCondSurvProb))
              }


res_true = foreach(pid=testData$testDs.id$P_ID, 
                   .packages = c("splines", "JMbayes"), .combine = 'rbind') %dopar%{
                     source("src/mdp/common/prediction.R")
                     
                     age = testData$testDs.id$Age[testData$testDs.id$P_ID == pid]
                     progression_time = testData$testDs.id$progression_time[testData$testDs.id$P_ID == pid]
                     true_b = unlist(testData$testDs.id[testData$testDs.id$P_ID == pid, 5:11])
                     
                     decision_epoch = 10
                     
                     survTimes = 1:10
                     survProbMatrix = getTrueFutureOutcomes(mvJoint_dre_psa_dre_value, age, 
                                                            true_b,0, Inf, survTimes, 
                                                            M = N_MCMC_ITER)$predicted_surv_prob
                     
                     meanSurvProb = rowMeans(survProbMatrix)
                     lowSurvProb = apply(survProbMatrix, 1, quantile, probs=0.025)
                     uppSurvProb = apply(survProbMatrix, 1, quantile, probs=0.975)
                     
                     wtPoints = getGaussianQuadWeightsPoints(c(0,10))
                     points = sort(wtPoints$points, decreasing = F)
                     wt = wtPoints$weights[order(wtPoints$points, decreasing = F)]
                     
                     survProbMatrix = getTrueFutureOutcomes(mvJoint_dre_psa_dre_value, age, 
                                                            true_b, 0, Inf, points, 
                                                            M = N_MCMC_ITER)$predicted_surv_prob
                     auc = sum(wt*rowMeans(survProbMatrix))
                     
                     return(data.frame('P_ID'=pid,'progression_time'=progression_time,
                                       'auc'=auc, 'survTimes'=survTimes, 'meanSurvProb'=meanSurvProb,
                                       'lowSurvProb'=lowSurvProb, 'uppSurvProb'=uppSurvProb))
                   }

stopCluster(ct)

res$correct = ifelse(res$progression_time <= 2, BIOPSY, WAIT)
res$correct_bin = ifelse(res$progression_time <= 2, 1, 0)

res_true$correct = ifelse(res_true$progression_time <= 2, BIOPSY, WAIT)
res_true$correct_bin = ifelse(res_true$progression_time <= 2, 1, 0)

combo = paste(res$correct, res$survTimes, sep="-")
by(data = res$meanSurvProb, INDICES = combo, FUN = summary)

res$decision_10perc = unlist(by(INDICES = res$P_ID, data = res, FUN = function(x){
  
  correct = ifelse(x$progression_time <=2, yes = BIOPSY, no=WAIT)
  done = ifelse(x$meanSurvProb[x$survTimes==2] <= 0.9, yes = BIOPSY, no = WAIT)
  
  return(paste0(correct, done))
}))

table(res$decision_10perc[!duplicated(res$P_ID)])

#Composite decisions
res$composite_dec_marginal = unlist(by(INDICES = res$P_ID, data = res,
                            FUN = function(x){
                              correct = ifelse(x$progression_time <=2, yes = BIOPSY, no=WAIT)
                              done = ifelse(x$meanSurvProb[x$survTimes==2] <= 0.9 &
                                              x$meanSurvProb[x$survTimes==5] <= 0.805, 
                                            yes = BIOPSY, no = WAIT)
                              
                              return(paste0(correct, done))
                            }))
table(res$composite_dec_marginal[!duplicated(res$P_ID)])

res$composite_dec_conditional = unlist(by(INDICES = res$P_ID, data = res,
                              FUN = function(x){
                                correct = ifelse(x$progression_time <=2, yes = BIOPSY, no=WAIT)
                                done = ifelse(x$meanSurvProb[x$survTimes==2] <= 0.9 &
                                                x$meanCondSurvProb[x$survTimes==5] <= 0.91, 
                                              yes = BIOPSY, no = WAIT)
                                return(paste0(correct, done))
                              }))
table(res$composite_dec_conditional[!duplicated(res$P_ID)])

res_wide = reshape(res, idvar = "P_ID", timevar = "survTimes", 
                   direction = "wide", v.names = c("meanSurvProb", "meanCondSurvProb"))

composite_benefit_rule = res_wide$P_ID %in% res_wide$P_ID[res_wide$decision_10perc != res_wide$composite_dec_conditional | 
                                                            res_wide$decision_10perc != res_wide$composite_dec_marginal]
View(res_wide[composite_benefit_rule,])
     
REWARDS = c(NA, NA, NA, NA)
names(REWARDS) = reward_names
REWARDS[TRUE_WAIT] = 0
REWARDS[FALSE_BIOPSY] = -2
REWARDS[TRUE_BIOPSY] = 2
REWARDS[FALSE_WAIT] = -7.090909

valueFuncs = function(meanSurvProb.2, meanSurvProb.5){
  V_WB = (1-meanSurvProb.2) * REWARDS[FALSE_WAIT] + meanSurvProb.2 * REWARDS[TRUE_WAIT] +
    (1-meanSurvProb.5) * REWARDS[TRUE_BIOPSY] + meanSurvProb.5 * REWARDS[FALSE_BIOPSY]
  
  V_WW = (1-meanSurvProb.2) * REWARDS[FALSE_WAIT] + meanSurvProb.2 * REWARDS[TRUE_WAIT] +
    (1-meanSurvProb.5) * REWARDS[FALSE_WAIT] + meanSurvProb.5 * REWARDS[TRUE_WAIT]
  
  V_BW = (1-meanSurvProb.2) * REWARDS[TRUE_BIOPSY] + meanSurvProb.2 * REWARDS[FALSE_BIOPSY] + 
    (meanSurvProb.2-meanSurvProb.5) * REWARDS[FALSE_WAIT] + meanSurvProb.5 * REWARDS[TRUE_WAIT]
  
  V_BB = (1-meanSurvProb.2) * REWARDS[TRUE_BIOPSY] + meanSurvProb.2 * REWARDS[FALSE_BIOPSY] + 
    (meanSurvProb.2-meanSurvProb.5) * REWARDS[TRUE_BIOPSY] + meanSurvProb.5 * REWARDS[FALSE_BIOPSY]
  
  return(max(c(V_BW, V_BW)) - max(c(V_WB, V_WW)))
}

res_wide[,c("VBW_minus_VWB")] = valueFuncs(res_wide$meanSurvProb.2, res_wide$meanSurvProb.5)

res_wide$mdp_reward_dec = ifelse(res_wide$VBW_minus_VWB > 0, yes = BIOPSY, no = WAIT)
res_wide$mdp_reward_dec = paste0(res_wide$correct, res_wide$mdp_reward_dec)

table(res_wide$mdp_reward_dec)
table(res_wide$decision_10perc)
table(res_wide$composite_dec_marginal)
table(res_wide$composite_dec_conditional)

###########################
# Trying to convert a composite decision into a MDP
###########################
plotDf = data.frame(expand.grid("TW"=0, "FB"=-2, "TB"=2, "FW"=seq(2, -30, length.out = 5000), 
                     "meanSurvProb.2"=c(0.95, 0.89, 0.85), 
                     "meanSurvProb.5"=c(0.85, 0.80, 0.75)))

plotDf[, c("VBW_minus_VWB")] = sapply(1:nrow(plotDf), function(i){
  REWARDS[TRUE_WAIT] <<- plotDf[i, TRUE_WAIT]
  REWARDS[FALSE_BIOPSY] <<- plotDf[i, FALSE_BIOPSY]
  REWARDS[TRUE_BIOPSY] <<- plotDf[i, TRUE_BIOPSY]
  REWARDS[FALSE_WAIT] <<- plotDf[i, FALSE_WAIT]
  valueFuncs(plotDf$meanSurvProb.2[i], plotDf$meanSurvProb.5[i])
})
plotDf$decision=ifelse(plotDf$VBW_minus_VWB > 0, BIOPSY, WAIT)

tt = expand.grid(seq(0,1, 0.01), seq(0,1, 0.01))
tt[, 3] = sapply(1:nrow(tt), function(i){
  valueFuncs(tt[i,1], tt[i,2])
})

uniroot(function(fw){
  REWARDS[FALSE_WAIT] <<- fw
  valueFuncs(0.9, 0.805)
}, c(0, -500))

ggplot(plotDf) + geom_line(aes(x=FW, y=VBW_minus_VWB, color=decision)) + 
  scale_color_manual(name="Decision",values =  c("red", "forestgreen")) + 
  facet_grid(paste("Two", meanSurvProb.2)~meanSurvProb.5) + theme_bw() + 
  theme(legend.position = "bottom")

############################
# Checking if survival curves are parallel for subjects whose year 2 survival prob is same
############################
res_true_2 = res_true[res_true$survTimes==2,]
res_true_2_sorted = res_true_2[order(res_true_2$meanSurvProb, decreasing = F),]
closest_pairs_id_1 = res_true_2_sorted$P_ID[order(diff(res_true_2_sorted$meanSurvProb), decreasing = F)]
closest_pairs_id_2 = res_true_2_sorted$P_ID[order(diff(res_true_2_sorted$meanSurvProb), decreasing = F)+1]

#Now I compare their entire curves
View(res_true[res_true$P_ID %in% c(closest_pairs_id_1[3], closest_pairs_id_2[3]),])
ggplot(data=res_true[res_true$P_ID %in% c(closest_pairs_id_1[3], closest_pairs_id_2[3]),]) +
  geom_line(aes(x=survTimes, y=meanSurvProb, 
                group=factor(P_ID), color=factor(P_ID))) +
  geom_ribbon(aes(x=survTimes, ymin=lowSurvProb, ymax=uppSurvProb,
                  group=factor(P_ID), fill=factor(P_ID)), alpha=0.5) + 
  ylim(0,1) + xlim(0,10)

#Can we say that two patients who have same surv prob at year 2.
#For them we look at their surv prob at year 5
#For each pair we expect that survival prob at year 5 should be ordered as their prog time
res_true_5 = res_true[res_true$survTimes==5,]

temp_res = c()
for(i in 1:length(closest_pairs_id_1)){
  prog_time_1 = res_true_5$progression_time[res_true_5$P_ID == closest_pairs_id_1[i]]
  prog_time_2 = res_true_5$progression_time[res_true_5$P_ID == closest_pairs_id_2[i]]
  
  surv_time_1 = res_true_5$meanSurvProb[res_true_5$P_ID == closest_pairs_id_1[i]]
  surv_time_2 = res_true_5$meanSurvProb[res_true_5$P_ID == closest_pairs_id_2[i]]
  
  temp_res[i] = paste(prog_time_1<=2, prog_time_2<=2,
                      surv_time_1 <= 0.9, surv_time_2 <=0.9,
                      abs(round(surv_time_1 - surv_time_2,1)), sep=" ; ")
}
table(temp_res)

#Now I apply the 10% rule and find patients who had an incorrect decision
res_true_2$correct = res_true_2$progression_time<=2
res_true_2$actual = res_true_2$meanSurvProb<0.9

table(paste(res_true_2$correct, res_true_2$actual))
incorrect_pid_2 = res_true_2$P_ID[res_true_2$correct != res_true_2$actual]

res_true_2_incorrect = res_true_2[res_true_2$P_ID %in% incorrect_pid_2,]
