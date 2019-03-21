# library(JMbayes)
# library(survival)
# library(splines)
# library(ggplot2)
# 
# #source the common methods for all algorithms
# source("src/mdp/common/simCommon.R")
# source("src/mdp/common/prediction.R")
# source("src/mdp/tree_search/forward_search_no_Y.R")
# 
# fitted_JM = jointModelData$mvJoint_dre_psa_simDs
# 
# patient_list = jointModelData$testData$testDs[jointModelData$testData$testDs$visitTimeYears<=2, ]
# patient_list = split(patient_list, patient_list$P_ID)
# 
# plist = c(1:250)[jointModelData$testData$testDs.id$progression_time>2]
# #plist = c(21,24, 28,3, 16, 1, 26, 45, 46, 35, 42, 4, 41)
# true_dec = ifelse(jointModelData$testData$testDs.id$progression_time<=2, "B","W")
# surv_dec = rep(NA, 250)
# my_dec = rep(NA, 250)
# for(i in plist){
#   set.seed(100)
#   times = seq(0, 10, by = 0.25)[-1]
#   fut = getExpectedFutureOutcomes(fitted_JM, patient_list[[i]], 0, Inf, 
#                                   survival_predict_times = times,
#                                   psa_predict_times = times,
#                                   dre_predict_times = times, M = 0)
#   meansurv = rowMeans(fut$predicted_surv_prob)
#   meanPSA = meanPSASlope = rowMeans(fut$predicted_psa)
#   meanPSASlope = rowMeans(fut$predicted_psa_slope)
#   meanPSAAccel = rowMeans(fut$predicted_psa_acceleration)
#   meanDRE = rowMeans(fut$predicted_dre_prob)
#   
#   p1 = ggplot() + geom_line(aes(x=times, y=meansurv)) + ylim(0,1) +ylab("surv")+
#     ggtitle(paste(i,"-", round(patient_list[[i]]$progression_time[1],1)))
#   p1_1 = ggplot() + geom_line(aes(x=times, y=meanPSA)) + ylab("value")+
#     ggtitle(paste(i,"-", round(patient_list[[i]]$progression_time[1],1)))
#   p2 = ggplot() + geom_line(aes(x=times, y=meanPSASlope)) + ylim(-0.5, 0.5) + ylab("Velo")+
#     ggtitle(paste(i,"-", round(patient_list[[i]]$progression_time[1],1)))
#   p3 = ggplot() + geom_line(aes(x=times, y=meanPSAAccel)) + ylab("accel") +
#     ggtitle(paste(i,"-", round(patient_list[[i]]$progression_time[1],1)))
#   p4 = ggplot() + geom_line(aes(x=times, y=meanDRE)) + ylim(0,1) + ylab("dreprob") +
#     ggtitle(paste(i,"-", round(patient_list[[i]]$progression_time[1],1)))
#   
#   # p1 = ggplot() + geom_line(aes(x=times, y=meansurv)) + ylim(0,1) + xlim(2, 10) +
#   #   geom_hline(yintercept = 0.9, color='red')
#   # p2 = ggplot() + geom_line(aes(x=times, y=meanPSA)) + ylim(-0.5, 0.5) + xlim(2,10)
#   # p3 = ggplot() + geom_line(aes(x=times, y=meanDRE)) + ylim(0,1)  + xlim(2,10)
#   p = ggpubr::ggarrange(p1,p1_1,p2,p3,p4, ncol=3, nrow=2)
#   print(p)
#   dec = readline()
#   if(dec=="e"){
#     break
#   }
#   #else{
#   #  surv_dec[i] = ifelse(meansurv[8]<=0.9, "B","W")
#   #  my_dec[i] = ifelse(meansurv[8]<=0.9 && (meanDRE[8]>=0.10 | meanPSA[8]>0.15), "B", "W")
#   #}
# }



#Levels
#1. reward row
#2. discount factors
#3. max biopsies: 1,2,5, Inf
#4. max depths: 0, 5 
discount_factors = as.factor(seq(0.5, 5, 0.5))
max_biopsies_vec = as.factor(c(1, 2, 5, Inf))
max_depths = as.factor(c(0,5))

reward_names = c(TRUE_BIOPSY, FALSE_BIOPSY, TRUE_WAIT, FALSE_WAIT)
# reward_matrix = as.matrix(expand.grid(seq(5,100,by=10),
#                                       -1, 1,
#                                       seq(-5,-100,by=-10)))
reward_matrix = as.matrix(expand.grid(5*(10^(-6:6)),
                                      -1, 1,
                                      -5*(10^(-6:6))))

colnames(reward_matrix) = reward_names

temp = vector("list", nrow(reward_matrix)*
                length(max_depths)*
                length(max_biopsies_vec)*
                length(discount_factors))

m = 1
for(i in 1:nrow(reward_matrix)){
  REWARD = reward_matrix[i,]
  
  for(j in 1:length(discount_factors)){
    for(k in 1:length(max_biopsies_vec)){
      for(l in 1:length(max_depths)){
        max_depth = max_depths[l]
        max_depth_number = as.numeric(as.character(max_depth))
        
        ds = sim_res[[i]][[j]][[k]][[max_depth_number + 1]]
        
        ds$max_depth = max_depth
        ds$max_biopsies = max_biopsies_vec[k]
        ds$discount_factor = discount_factors[j]
        ds$R_TB = as.factor(reward_matrix[i, TRUE_BIOPSY])
        ds$R_TW = reward_matrix[i, TRUE_WAIT]
        ds$R_FB = reward_matrix[i, FALSE_BIOPSY]
        ds$R_FW = as.factor(reward_matrix[i, FALSE_WAIT])
        
        temp[[m]] = ds
        m = m+1
      }
    }
  }
}

combined_sim_res = do.call('rbind', temp)
combined_sim_res$nb = as.numeric(combined_sim_res$nb)

effects = glmer(nb~max_depth + max_biopsies + (1|P_ID), data=combined_sim_res,
      family = binomial)

View(combined_sim_res[combined_sim_res$max_biopsies==1,])
