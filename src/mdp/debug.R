res = list.files("Rdata/mdp/DF_pt9/", full.names = T)

finalRes = do.call(rbind, lapply(res, FUN = function(filePath){
  threshold = substr(filePath, 26,30)
  
  load(filePath)
  
  decision_epochs = c(1,1.25,1.5,1.75,2,2.5,3,3.5,4,4.5,5)
  resMatrix = do.call(rbind, lapply(as.character(decision_epochs), 
                      function(max_decision_epoch_str){
                        do.call('c',by(dataset.id$`1`[,max_decision_epoch_str], 
                                       INDICES = dataset.id$`1`$progression_time<=1, table))
                      }))
  
  return(data.frame(threshold=threshold, 
                    decision_epochs, 
                    nBB=resMatrix[,1], 
                    nGW=resMatrix[,2],
                    nGB=resMatrix[,3], 
                    nBW=resMatrix[,4]))
  
}))

lapply(sim_res[[1]]$`2`, function(x){
  action_taken = sapply(x, function(y){y$optimal_action})
  optimal_action = ifelse(jointModelData$testData$testDs.id$progression_time<=1, "B", "W")
  return(c("nTB"=sum(action_taken=="B" & optimal_action=="B"),
           "nFB"=sum(action_taken=="B" & optimal_action=="W"),
           "nTW"=sum(action_taken=="W" & optimal_action=="W"),
           "nFW"=sum(action_taken=="W" & optimal_action=="B")))
})
