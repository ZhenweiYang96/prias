library(JMbayes)
library(survival)
library(splines)
library(ggplot2)

#source the common methods for all algorithms
source("src/mdp/common/simCommon.R")
source("src/mdp/common/prediction.R")

max_cores = 3
N_MCMC_ITER = 300

file = list.files(path='C://Users//838035//prias//Rdata//mdp//final_res', full.names = T)[1]
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

          return(data.frame('P_ID'=pid,'progression_time'=progression_time,
                            'auc'=auc, 'survTimes'=survTimes, 'meanSurvProb'=meanSurvProb,
                   'lowSurvProb'=lowSurvProb, 'uppSurvProb'=uppSurvProb))
        }


res_true = foreach(pid=testData$testDs.id$P_ID, 
              .packages = c("splines", "JMbayes"), .combine = 'rbind') %dopar%{
                source("src/mdp/common/prediction.R")
                
                patient_df = testData$testDs[testData$testDs$P_ID==pid,]
                progression_time = patient_df$progression_time[1]
                true_b = unlist(testData$testDs.id[testData$testDs.id$P_ID == pid, 5:11])
                
                decision_epoch = 10
                pat_data = patient_df[patient_df$visitTimeYears <= decision_epoch,]
                
                survTimes = 1:10
                survProbMatrix = getTrueFutureOutcomes(mvJoint_dre_psa_dre_value, pat_data, 
                                                       true_b,0, Inf, survTimes, 
                                                       M = N_MCMC_ITER)$predicted_surv_prob
                
                meanSurvProb = rowMeans(survProbMatrix)
                lowSurvProb = apply(survProbMatrix, 1, quantile, probs=0.025)
                uppSurvProb = apply(survProbMatrix, 1, quantile, probs=0.975)
                
                wtPoints = getGaussianQuadWeightsPoints(c(0,10))
                points = sort(wtPoints$points, decreasing = F)
                wt = wtPoints$weights[order(wtPoints$points, decreasing = F)]
                
                survProbMatrix = getTrueFutureOutcomes(mvJoint_dre_psa_dre_value, pat_data, 
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

res_true$correct = ifelse(res$progression_time <= 2, BIOPSY, WAIT)
res_true$correct_bin = ifelse(res$progression_time <= 2, 1, 0)


res_wide = reshape(res, idvar = "P_ID", timevar = "survTimes", direction = "wide",v.names = c("meanSurvProb"))

combo = paste(res$correct, res$survTimes, sep="-")
by(data = res$meanSurvProb, INDICES = combo, FUN = summary)


perf = res_true[!duplicated(res_true$P_ID),]

perf$decision_10perc_true = by(INDICES = res_true$P_ID, data = res_true, FUN = function(x){
  
  correct = ifelse(x$progression_time[1] <=2, yes = BIOPSY, no=WAIT)
  done = ifelse(x$meanSurvProb[x$survTimes==2] <= 0.9, yes = BIOPSY, no = WAIT)
  
  return(paste0(correct, done))
})

perf$decision_10perc = by(INDICES = res$P_ID, data = res, FUN = function(x){
  
  correct = ifelse(x$progression_time[1] <=2, yes = BIOPSY, no=WAIT)
  done = ifelse(x$meanSurvProb[x$survTimes==2] <= 0.9, yes = BIOPSY, no = WAIT)
  
  return(paste0(correct, done))
})


res$decision_10perc = unlist(by(INDICES = res$P_ID, data = res, FUN = function(x){
  
  correct = ifelse(x$progression_time <=2, yes = BIOPSY, no=WAIT)
  done = ifelse(x$meanSurvProb[x$survTimes==2] <= 0.9, yes = BIOPSY, no = WAIT)
  
  return(paste0(correct, done))
}))

res_BB = res[res$P_ID %in% perf$P_ID[perf$decision_10perc=="BB"],]
res_WB = res[res$P_ID %in% perf$P_ID[perf$decision_10perc=="WB"],]
res_WW = res[res$P_ID %in% perf$P_ID[perf$decision_10perc=="WW"],]
res_BW = res[res$P_ID %in% perf$P_ID[perf$decision_10perc=="BW"],]

perf$decision_10perc_and_315perc = by(INDICES = res$P_ID, data = res, FUN = function(x){
  correct = ifelse(x$progression_time[1] <=2, yes = BIOPSY, no=WAIT)
  done = ifelse(x$meanSurvProb[x$survTimes==2] <= 0.9 &
                  x$meanSurvProb[x$survTimes==5] <= 0.7, yes = BIOPSY, no = WAIT)
  
  return(paste0(correct, done))
})
