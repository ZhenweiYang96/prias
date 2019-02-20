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
