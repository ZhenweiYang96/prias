load("/home/a_tomer/Data/mdp/schedule_by_reward/sim_res_10_36.Rdata")

max_depths = 0:4

reward_names = c(TRUE_BIOPSY, FALSE_BIOPSY, TRUE_WAIT, FALSE_WAIT)
reward_matrix = as.matrix(expand.grid(seq(5,100,by=10),
                                      -1, 1, 
                                      seq(-5,-100,by=-10)))
colnames(reward_matrix) = reward_names

orig_res = jointModelData$testData$scheduleResults

perf_10_perc = orig_res[orig_res$methodName %in% c("Risk (10%)","Risk (15%)"),]

perf_10_perc$nb[perf_10_perc$offset < 0] = perf_10_perc$nb[perf_10_perc$offset < 0] + 1
perf_10_perc$offset[perf_10_perc$offset < 0] = 10 - perf_10_perc$progression_time[perf_10_perc$offset < 0]


sim_res = lapply(sim_res, function(x){
  lapply(x, function(y){
    y$nb[y$offset < 0] = y$nb[y$offset < 0] + 1
    y$offset[y$offset < 0] = 10 - y$progression_time[y$offset < 0]
    return(y)
  })
})

for(i in temp){
  #for(i in 1){
  final_res = perf_10_perc
  final_res = cbind(final_res[,-4],final_res[,c("methodName")])
  colnames(final_res)[6] = "methodName" 
  
  for(max_depth in max_depths){
    #for(max_depth in 0){
    new_data = sim_res[[i]][[max_depth + 1]]
    new_data$methodName = paste0("Rwd", "_", reward_matrix[i,1], 
                                 "_", reward_matrix[i,4], "_", max_depth)
    final_res = rbind(final_res, new_data)
  }
  
  a=ggplot(data=final_res[final_res$progression_time==10,]) + 
    geom_boxplot(aes(y=nb, x=methodName, fill=grepl("Risk", methodName))) + 
    ggtitle(paste("NB-NP", reward_matrix[i,"TB"], reward_matrix[i,"FW"])) +
    theme(legend.position = "none")
  
  b = ggplot(data=final_res[final_res$progression_time<10,]) + 
    geom_boxplot(aes(y=nb, x=methodName, fill=grepl("Risk", methodName))) + 
    ggtitle(paste("NB-P", reward_matrix[i,"TB"], reward_matrix[i,"FW"])) +
    theme(legend.position = "none")
  
  c=ggplot(data=final_res[final_res$progression_time<10,]) + 
    geom_boxplot(aes(y=offset, x=methodName, fill=grepl("Risk", methodName))) + 
    ggtitle(paste("Offset-P", reward_matrix[i,"TB"], reward_matrix[i,"FW"])) +
    theme(legend.position = "none")
  
  now = ggpubr::ggarrange(a,b,c, ncol = 3, nrow=1)
   if(i > 90){
     print(ggpubr::ggarrange(prev, now, ncol=1, nrow = 2))
   }

   prev = now
  #print(now)
}


##########
#2 D plots
##########

for(max_depth in 4){
  size = sapply(sim_res, function(x){
    median(x[[max_depth + 1]]$nb)
  })
}

ggplot() + geom_point(aes(x=reward_matrix[,1], y=reward_matrix[,4], color=size))
