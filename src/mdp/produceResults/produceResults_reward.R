load("/home/a_tomer/Data/mdp/schedule_by_reward/sim_res_10_36.Rdata")

max_depths = 0:4

reward_names = c(TRUE_BIOPSY, FALSE_BIOPSY, TRUE_WAIT, FALSE_WAIT)
# reward_matrix = as.matrix(expand.grid(c(5,10,15,20,25,30),
#                                       -1, 1, 
#                                       c(-5,-10,-15,-20,-25,-30)))

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

#for(i in 1:length(sim_res)){
for(i in 1){
  
  final_res = perf_10_perc
  final_res = cbind(final_res[,-4],final_res[,c("methodName")])
  colnames(final_res)[6] = "methodName" 
  
  for(max_depth in 0:2){
  #for(max_depth in 0){
    new_data = sim_res[[i]][[max_depth + 1]]
    new_data$methodName = paste0("Rwd", "_", reward_matrix[i,1], 
                                 "_", reward_matrix[i,4], "_", max_depth)
    final_res = rbind(final_res, new_data)
  }
  
  a=ggplot(data=final_res[final_res$progression_time==10,]) + geom_boxplot(aes(y=nb, x=methodName)) + ggtitle(paste("NB-NP", i))
  b = ggplot(data=final_res[final_res$progression_time<10,]) + geom_boxplot(aes(y=nb, x=methodName)) + ggtitle(paste("NB", i))
  c=ggplot(data=final_res[final_res$progression_time<10,]) + geom_boxplot(aes(y=offset, x=methodName)) + ggtitle(paste("Offset", i))
  print(ggpubr::ggarrange(a,b,c, ncol = 3, nrow=1))
}


