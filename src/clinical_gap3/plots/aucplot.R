center_names = sapply(list.files("Rdata/gap3/PRIAS_2019/auc_only_psa/"), function(x){
  sub("*.Rdata", "", sub("auc_list_*", "", x))
})

auc_file_names = list.files("Rdata/gap3/PRIAS_2019/auc_only_psa/", full.names = T)
prederr_file_names = list.files("Rdata/gap3/PRIAS_2019/prederr_only_psa/", full.names = T)

auc_onlypsa_file_names = list.files("Rdata/gap3/PRIAS_2016/auc_only_psa/", full.names = T)
prederr_onlypsa_file_names = list.files("Rdata/gap3/PRIAS_2016/prederr_only_psa/", full.names = T)

plotDf = do.call('rbind', lapply(1:length(center_names), function(i){
  print(i)
  
  load(auc_file_names[i])
  aucs = sapply(auc_list, function(x){x$auc})
  auc_nrs = sapply(auc_list, function(x){
    if(is.null(x$nr)){
      NA
    }else{
      x$nr
    }
  })
  
  load(auc_onlypsa_file_names[i])
  aucs_onlypsa = sapply(auc_list, function(x){x$auc})
  auc_nrs_onlypsa = sapply(auc_list, function(x){
    if(is.null(x$nr)){
      NA
    }else{
      x$nr
    }
  })
  
  load(prederr_file_names[i])
  prederrs = sqrt(sapply(prederr_list, function(x){x$prederr}))
  prederr_nrs = sapply(prederr_list, function(x){x$nr})
  
  load(prederr_onlypsa_file_names[i])
  prederrs_onlypsa = sqrt(sapply(prederr_list, function(x){x$prederr}))
  prederr_nrs_onlypsa = sapply(prederr_list, function(x){x$nr})
  
  data.frame('Center'=center_names[i], 
             'auc'=aucs, 'auc_nr'=auc_nrs,
             'auc_onlypsa'=aucs_onlypsa, 'auc_nr_onlypsa'=auc_nrs_onlypsa,
             'prederr'=prederrs, 'prederr_nr'=prederr_nrs,
             'prederr_onlypsa'=prederrs_onlypsa, 'prederr_nr'=prederr_nrs_onlypsa,
             'time'=1:10)
}))

pdf("Rdata/gap3/PRIAS_2019/plots.pdf", width = 12, height = 7)

lapply(center_names[!(center_names %in% c("prias_orig", "PRIAS", "Gothenburg"))], function(x){
  temp = plotDf[plotDf$Center %in% c("PRIAS", x),]
  temp$Center = as.character(temp$Center)
  
  temp$Center[temp$Center == "PRIAS"] = "1. PRIAS"
  temp$Center[temp$Center == x] = paste0("2. ", x)
  
  auc = ggplot(data = temp) +
    geom_line(aes(x=time, y=auc, color=Center)) +
    #geom_line(aes(x=time, y=auc_onlypsa, color=Center), linetype='dashed') + 
    ylab("AUC") + xlab("Time (years)")+
    scale_x_continuous(breaks = c(1:10)) +
    theme_bw()+
    theme(text=element_text(size=15)) + 
    ylim(0,1) + geom_hline(yintercept = 0.5) + 
    geom_vline(xintercept = 5)
  #geom_label(aes(x=time, y=auc, label=auc_nr)) 
  
  prederr = ggplot(data = temp) + 
    geom_line(aes(x=time, y=prederr, color=Center)) +
    #geom_line(aes(x=time, y=prederr_onlypsa, color=Center), linetype='dashed') + 
    scale_x_continuous(breaks = c(1:10)) + 
    ylab("Sqrt (Prediction Error or MSE)") + xlab("Time (years)")+
    theme_bw()+
    theme(text=element_text(size=15)) + 
    geom_vline(xintercept = 5) + 
    ylim(0,1) 
   #geom_label(aes(x=time, y=auc, label=prederr_nr))
  
  ggpubr::ggarrange(auc, prederr, ncol = 2, common.legend = T)
})

dev.off()
