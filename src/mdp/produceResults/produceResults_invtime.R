load("/home/a_tomer/Data/mdp/schedule_by_invtime/sim_res_10_risk15pc.Rdata")

max_depths = 0:4

reward_names = c(TRUE_BIOPSY, FALSE_BIOPSY, TRUE_WAIT, FALSE_WAIT)
reward_matrix = as.matrix(expand.grid(c(5,10,15,20,25,30),
                                      -1, 1, 
                                      c(-5,-10,-15,-20,-25,-30)))
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

final_res = perf_10_perc
final_res = cbind(final_res[,-4],final_res[,c("methodName")])
colnames(final_res)[6] = "methodName" 
for(i in 1:length(sim_res)){
#for(i in 1:1){
  for(max_depth in 0){
  #for(max_depth in 0){
    new_data = sim_res[[i]][[max_depth + 1]]
    new_data$methodName = paste0("Rwd", "_", names(sim_res)[i], "_", max_depth)
    final_res = rbind(final_res, new_data)
  }
}

a=ggplot(data=final_res) + geom_boxplot(aes(y=nb, x=methodName))
b=ggplot(data=final_res) + geom_boxplot(aes(y=offset, x=methodName))
ggpubr::ggarrange(a,b)
