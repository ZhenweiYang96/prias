combined_list = list()
k = 1

for(t in 1:length(thresholds)){
  threshold = thresholds[t]
  
  #for(df in 1:length(discount_factors)){
  for(df in 1){
    DISCOUNT_FACTOR = discount_factors[df]
    
    for(v in 1:length(value_functions)){
      value = value_functions[v]
      
      for(bs in 1:length(slopes)){
        biopsy_slope = slopes[bs]
        
        biopsy_intercept = value - biopsy_slope * threshold
        
        for(ws in 1:length(slopes)){
          wait_slope = slopes[ws]
          
          if(wait_slope != biopsy_slope){
            wait_intercept = value - wait_slope*threshold
            
            REWARDS = thresholdToReward(NA, biopsy_intercept, biopsy_slope, wait_intercept, wait_slope)
            
            combined_list[[k]] = sim_res[[as.character(threshold)
                                          ]][[as.character(DISCOUNT_FACTOR)
                                              ]][[as.character(value)
                                                  ]][[as.character(biopsy_slope)
                                                      ]][[as.character(wait_slope)]]
            
            combined_list[[k]]$correct_action = ifelse(combined_list[[k]]$progression_time<=2, BIOPSY, WAIT)
            combined_list[[k]]$wait_intercept = factor(wait_intercept)
            combined_list[[k]]$wait_slope = factor(wait_slope)
            combined_list[[k]]$biopsy_intercept = factor(biopsy_intercept)
            combined_list[[k]]$biopsy_slope = factor(biopsy_slope)
            combined_list[[k]]$value = factor(value)
            combined_list[[k]]$DISCOUNT_FACTOR = factor(DISCOUNT_FACTOR)
            combined_list[[k]]$threshold = factor(threshold)
            combined_list[[k]][, c("R_TB", "R_FB", "R_TW", "R_FW")] = rep(REWARDS, each=100)
              
            k = k+1
          }
        }
      }
    }
  }
}

combined_res = do.call('rbind', combined_list)
combined_summary = do.call('rbind', lapply(combined_list, function(x){
  ret = x[1,-c(1:6)]
  ret[,c("TBR", "FBR", "TWR", "FWR")] = c(sum(x$optimal_action==BIOPSY & x$correct_action==BIOPSY),
                                              sum(x$optimal_action==BIOPSY & x$correct_action==WAIT),
                                              sum(x$optimal_action==WAIT & x$correct_action==WAIT),
                                              sum(x$optimal_action==WAIT & x$correct_action==BIOPSY))/nrow(x)
  return(ret)
}))

ggplot(data=combined_summary[combined_summary$threshold==0.1,]) + 
  geom_point(aes(x=FBR/(FBR+TWR),y=TBR/(TBR+FWR), shape=DISCOUNT_FACTOR, color=factor(slope_diff)), size=4, alpha=1) + 
  xlab("FBR") + ylab("TBR") + xlim(0,1) + ylim(0,1) + 
  geom_abline(slope=1, intercept = 0, color='red') + 
  facet_wrap(~intercept_diff)  + theme(legend.position = "bottom") 
