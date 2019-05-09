final_res = list()

for(i in 1:nrow(reward_matrix)){
  print(paste0("Running for reward row (",i,"): ", paste(reward_matrix[i,], collapse = ' ')))
  
  REWARDS = reward_matrix[i,]
  
  for(j in 1:length(discount_factors)){
    DISCOUNT_FACTOR = discount_factors[j]
    print(paste("Running for discount factor:", DISCOUNT_FACTOR))
    
    for(q in 1:length(max_biopsies_vec)){
      max_biopsies = max_biopsies_vec[q]
      print(paste("Running for max_biopsies:", max_biopsies))
      
      for(max_depth in max_depths){
        
        print(paste("Max depth", max_depth))      
        pat_subset = 1:100
        
        res = sim_res[[i]][[j]][[q]][[max_depth + 1]]
        res$nb[res$offset<0] = res$nb[res$offset<0] + 1
        res$offset[res$offset<0] = 10 - res$progression_time[res$offset<0]
        
        meanNb = mean(res$nb)
        medianNb = median(res$nb)
        IQRNb = IQR(res$nb)
        medianNbProg = median(res$nb[res$progression_time < 10])
        IQRNbProg = IQR(res$nb[res$progression_time < 10])
        
        meanOffset = mean(res$offset)
        medianOffset = median(res$offset)
        IQROffset = IQR(res$offset)
        medianOffsetProg = median(res$offset[res$progression_time < 10])
        IQROffsetProg = IQR(res$offset[res$progression_time < 10])
        
        final_res[[length(final_res) + 1]] = c(meanNb=meanNb,medianNb=medianNb,
                                               IQRNb=IQRNb,medianNbProg=medianNbProg,
                                               IQRNbProg=IQRNbProg, meanOffset=meanOffset,
                                               medianOffset=medianOffset,IQROffset=IQROffset,
                                               medianOffsetProg=medianOffsetProg,
                                               IQROffsetProg=IQROffsetProg, 
                                               DISCOUNT_FACTOR=DISCOUNT_FACTOR, 
                                               max_depth=max_depth, REWARDS=REWARDS)
      }
    }
  }
}

combined_summary = data.frame(do.call('rbind', final_res))
