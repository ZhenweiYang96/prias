auc_files = list.files("Rdata/gap3/PRIAS_2019/auc_prederr/temp/tt/", pattern = "auc", full.names = T)
pred_err_files = list.files("Rdata/gap3/PRIAS_2019/auc_prederr/temp/tt/", pattern = "pred", full.names = T)

temp = vector("list", length = 25)
lapply(1:length(temp), FUN = function(i){
  load(auc_files[i])
  temp[[i]]$auc_list <<- auc_prederr_bs[[1]]$auc_list
  
  load(pred_err_files[i])
  temp[[i]]$prederr_list <<- auc_prederr_bs[[1]]$prederr_list
})

auc_prederr_bs = temp
