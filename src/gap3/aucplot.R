center_names = sapply(list.files("Rdata/gap3/auc/"), function(x){
  sub("*.Rdata", "", sub("auc_list_*", "", x))
})
auc_file_names = list.files("Rdata/gap3/auc/", full.names = T)
prederr_file_names = list.files("Rdata/gap3/prederr/", full.names = T)

plotDf = do.call('rbind', lapply(1:length(center_names), function(i){
  print(i)
  load(auc_file_names[i])
  load(prederr_file_names[i])
  
  aucs = sapply(auc_list, function(x){x$auc})
  auc_nrs = sapply(auc_list, function(x){
    if(is.null(x$nr)){
      NA
    }else{
      x$nr
    }
  })
  prederrs = sqrt(sapply(prederr_list, function(x){x$prederr}))
  prederr_nrs = sapply(prederr_list, function(x){x$nr})
  
  data.frame('center'=center_names[i], 
               'auc'=aucs, 'auc_nr'=auc_nrs,
             'prederr'=prederrs, 'prederr_nr'=prederr_nrs,
               'time'=1:10)
}))


lapply(center_names[center_names!=c("prias_orig", "PRIAS")], function(x){
  temp = plotDf[plotDf$center %in% c("prias_orig", x),]
  auc = ggplot(data = temp) + geom_line(aes(x=time, y=auc, color=center)) +
    ylim(0,1) + geom_hline(yintercept = 0.5)
    #geom_label(aes(x=time, y=auc, label=auc_nr)) +
  prederr = ggplot(data = temp) + geom_line(aes(x=time, y=prederr, color=center)) +
    ylim(0,1) + geom_hline(yintercept = 0.5)
  #+ geom_label(aes(x=time, y=auc, label=prederr_nr))
  
  print(ggpubr::ggarrange(auc, prederr, ncol = 2, common.legend = T))
})
