# auc_files = list.files("Rdata/gap3/PRIAS_2019/auc_prederr/temp/", pattern = "auc", full.names = T)
# pred_err_files = list.files("Rdata/gap3/PRIAS_2019/auc_prederr/temp/", pattern = "pred", full.names = T)
# 
# temp = vector("list", length = length(auc_files))
# lapply(1:length(temp), FUN = function(i){
#   load(auc_files[i])
#   temp[[i]]$auc_list <<- auc_prederr_bs[[1]]$auc_list
#   
#   load(pred_err_files[i])
#   temp[[i]]$prederr_list <<- auc_prederr_bs[[1]]$prederr_list
# })
# 
# auc_prederr_bs = temp
# sapply(temp, function(x){
#   sapply(x$auc_list, function(x){
#     x$auc
#   })
# })
# 
# sapply(temp, function(x){
#   sapply(x$prederr_list, function(x){
#     x$prederr
#   })
# })




fileNames_all = list.files("~/Data/AUC_Brier/all/", full.names = T)

auc_prederr_all = lapply(fileNames_all, function(x){
  print(x)
  load(x)
  bootstrapdata$mvJoint = NULL
  return(bootstrapdata)
})

fileNames_nodre = list.files("~/Data/AUC_Brier/no_dre/", full.names = T)
auc_prederr_nodre = lapply(fileNames_nodre, function(x){
  print(x)
  load(x)
  bootstrapdata$mvJoint = NULL
  return(bootstrapdata)
})

save.image(file="/home/a_tomer/Google Drive/PhD/src/prias/Rdata/decision_analytic/AUC_Brier/auc_summary.Rdata")

auc_prederr_all = lapply(auc_prederr_all, function(x){
  resMatrix = matrix(NA, ncol = 5, nrow=4)
  
  resMatrix[1,] = sapply(x$aucAll_bb, function(x){x$auc})
  resMatrix[2,] = sapply(x$aucAll_bo, function(x){x$auc})
  resMatrix[3,] = sapply(x$prederrAll_bb, function(x){x$prederr})
  resMatrix[4,] = sapply(x$prederrAll_bo, function(x){x$prederr})
  
  rownames(resMatrix) = c("aucAll_bb", "aucAll_bo", "prederrAll_bb", "prederrAll_bo")
  colnames(resMatrix) = 0:4
  
  return(list(seed = x$seed, resMatrix = resMatrix))
})

auc_prederr_nodre_nonull = sapply(auc_prederr_nodre, function(x){!is.null(x$auc_nodre_bb)})

auc_prederr_nodre = lapply(auc_prederr_nodre[auc_prederr_nodre_nonull], function(x){
  resMatrix = matrix(NA, ncol = 5, nrow=4)
  
  resMatrix[1,] = sapply(x$auc_nodre_bb, function(x){x$auc})
  resMatrix[2,] = sapply(x$auc_nodre_bo, function(x){x$auc})
  resMatrix[3,] = sapply(x$prederr_nodre_bb, function(x){x$prederr})
  resMatrix[4,] = sapply(x$prederr_nodre_bo, function(x){x$prederr})
  
  rownames(resMatrix) = c("auc_nodre_bb", "auc_nodre_bo", "prederr_nodre_bb", "prederr_nodre_bo")
  colnames(resMatrix) = 0:4
  
  return(list(seed = x$seed, resMatrix = resMatrix))
})

######################
matrices_all = lapply(auc_prederr_all, function(x){
  resMatrix = x$resMatrix
  
  resMatrix = rbind(resMatrix, x$resMatrix[1,]-x$resMatrix[2,])
  resMatrix = rbind(resMatrix, x$resMatrix[3,]-x$resMatrix[4,])
  rownames(resMatrix)[5:6] = c("O_auc", "O_prederr")
  return(resMatrix)
})

matrices_nodre = lapply(auc_prederr_nodre, function(x){
  resMatrix = x$resMatrix
  
  resMatrix = rbind(resMatrix, x$resMatrix[1,]-x$resMatrix[2,])
  resMatrix = rbind(resMatrix, x$resMatrix[3,]-x$resMatrix[4,])
  rownames(resMatrix)[5:6] = c("O_auc", "O_prederr")
  return(resMatrix)
})

matrix_all_auc_bb = t(sapply(matrices_all, function(x){
  x[1,]
}))

matrix_all_auc_bo = t(sapply(matrices_all, function(x){
  x[2,]
}))

matrix_all_prederr_bb = t(sapply(matrices_all, function(x){
  x[3,]
}))

matrix_all_prederr_bo = t(sapply(matrices_all, function(x){
  x[4,]
}))

##

matrix_nodre_auc_bb = t(sapply(matrices_nodre, function(x){
  x[1,]
}))

matrix_nodre_auc_bo = t(sapply(matrices_nodre, function(x){
  x[2,]
}))

matrix_nodre_prederr_bb = t(sapply(matrices_nodre, function(x){
  x[3,]
}))

matrix_nodre_prederr_bo = t(sapply(matrices_nodre, function(x){
  x[4,]
}))

auc_final_all_0_1 = auc_0_1$auc - sapply(matrices_all, function(x){x[5,1]})
auc_final_all_1_2 = auc_1_2$auc - sapply(matrices_all, function(x){x[5,2]})
auc_final_all_2_3 = auc_2_3$auc - sapply(matrices_all, function(x){x[5,3]})
auc_final_all_3_4 = auc_3_4$auc - sapply(matrices_all, function(x){x[5,4]})
auc_final_all_4_5 = auc_4_5$auc - sapply(matrices_all, function(x){x[5,5]})

prederr_final_all_0_1 = prederr_0_1$prederr - sapply(matrices_all, function(x){x[6,1]})
prederr_final_all_1_2 = prederr_1_2$prederr - sapply(matrices_all, function(x){x[6,2]})
prederr_final_all_2_3 = prederr_2_3$prederr - sapply(matrices_all, function(x){x[6,3]})
prederr_final_all_3_4 = prederr_3_4$prederr - sapply(matrices_all, function(x){x[6,4]})
prederr_final_all_4_5 = prederr_4_5$prederr - sapply(matrices_all, function(x){x[6,5]})


auc_final_nodre_0_1 = auc_0_1_nodre$auc - sapply(matrices_nodre, function(x){x[5,1]})
auc_final_nodre_1_2 = auc_1_2_nodre$auc - sapply(matrices_nodre, function(x){x[5,2]})
auc_final_nodre_2_3 = auc_2_3_nodre$auc - sapply(matrices_nodre, function(x){x[5,3]})
auc_final_nodre_3_4 = auc_3_4_nodre$auc - sapply(matrices_nodre, function(x){x[5,4]})
auc_final_nodre_4_5 = auc_4_5_nodre$auc - sapply(matrices_nodre, function(x){x[5,5]})

prederr_final_nodre_0_1 = prederr_0_1_nodre$prederr - sapply(matrices_nodre, function(x){x[6,1]})
prederr_final_nodre_1_2 = prederr_1_2_nodre$prederr - sapply(matrices_nodre, function(x){x[6,2]})
prederr_final_nodre_2_3 = prederr_2_3_nodre$prederr - sapply(matrices_nodre, function(x){x[6,3]})
prederr_final_nodre_3_4 = prederr_3_4_nodre$prederr - sapply(matrices_nodre, function(x){x[6,4]})
prederr_final_nodre_4_5 = prederr_4_5_nodre$prederr - sapply(matrices_nodre, function(x){x[6,5]})

rm(list = setdiff(ls(), ls()[grep("_final+", ls(), perl=TRUE)]))

############
auc_graph = ggplot() + 
  geom_line(aes(x=1:5, y=c(mean(auc_final_all_0_1),
                           mean(auc_final_all_1_2),
                           mean(auc_final_all_2_3),
                           mean(auc_final_all_3_4),
                           mean(auc_final_all_4_5)), color="Both PSA and DRE")) +
  geom_ribbon(aes(x=1:5, ymin=sapply(list(auc_final_all_0_1,
                            auc_final_all_1_2, auc_final_all_2_3,
                            auc_final_all_3_4, auc_final_all_4_5), quantile, probs=0.025),
                  ymax = sapply(list(auc_final_all_0_1,
                              auc_final_all_1_2, auc_final_all_2_3,
                              auc_final_all_3_4, auc_final_all_4_5), quantile, probs=0.975)),
              fill='firebrick1', alpha=0.2) +
  geom_line(aes(x=1:5, y=c(mean(auc_final_nodre_0_1),
                           mean(auc_final_nodre_1_2),
                           mean(auc_final_nodre_2_3),
                           mean(auc_final_nodre_3_4),
                           mean(auc_final_nodre_4_5)), color="Only PSA")) +
  geom_ribbon(aes(x=1:5, ymin=sapply(list(auc_final_nodre_0_1,
                                          auc_final_nodre_1_2, auc_final_nodre_2_3,
                                          auc_final_nodre_3_4, auc_final_nodre_4_5), quantile, probs=0.025),
                  ymax = sapply(list(auc_final_nodre_0_1,
                                     auc_final_nodre_1_2, auc_final_nodre_2_3,
                                     auc_final_nodre_3_4, auc_final_nodre_4_5), quantile, probs=0.975)),
              fill='lightblue', alpha=0.5) +
  scale_color_manual(name="", values = c('red3', 'dodgerblue4')) + 
  xlab('Follow up time (years)') + ylab("AUC") + ylim(0,1) +
  theme_bw() + theme(legend.title = element_blank(), text=element_text(size=15),
                     legend.position = "bottom", legend.direction = "horizontal")

prederr_graph = ggplot() + 
  geom_line(aes(x=1:5, y=c(mean(prederr_final_all_0_1),
                           mean(prederr_final_all_1_2),
                           mean(prederr_final_all_2_3),
                           mean(prederr_final_all_3_4),
                           mean(prederr_final_all_4_5)), color="Both PSA and DRE")) +
  geom_ribbon(aes(x=1:5, ymin=sapply(list(prederr_final_all_0_1,
                                          prederr_final_all_1_2, prederr_final_all_2_3,
                                          prederr_final_all_3_4, prederr_final_all_4_5), quantile, probs=0.025),
                  ymax = sapply(list(prederr_final_all_0_1,
                                     prederr_final_all_1_2, prederr_final_all_2_3,
                                     prederr_final_all_3_4, prederr_final_all_4_5), quantile, probs=0.975)),
              fill='firebrick1', alpha=0.2) +
  geom_line(aes(x=1:5, y=c(mean(prederr_final_nodre_0_1),
                           mean(prederr_final_nodre_1_2),
                           mean(prederr_final_nodre_2_3),
                           mean(prederr_final_nodre_3_4),
                           mean(prederr_final_nodre_4_5)), color="Only PSA")) +
  geom_ribbon(aes(x=1:5, ymin=sapply(list(prederr_final_nodre_0_1,
                                          prederr_final_nodre_1_2, prederr_final_nodre_2_3,
                                          prederr_final_nodre_3_4, prederr_final_nodre_4_5), quantile, probs=0.025),
                  ymax = sapply(list(prederr_final_nodre_0_1,
                                     prederr_final_nodre_1_2, prederr_final_nodre_2_3,
                                     prederr_final_nodre_3_4, prederr_final_nodre_4_5), quantile, probs=0.975)),
              fill='lightblue', alpha=0.5) +
  scale_color_manual(name="", values = c('red3', 'dodgerblue4')) + 
  xlab('Follow up time (years)') + ylab("Prediction Error") + ylim(0, 0.16) +
  theme_bw() + theme(legend.title = element_blank(), text=element_text(size=15),
                     legend.position = "bottom", legend.direction = "horizontal")
################
ggplot() +
  geom_point(aes(x=0:4,y=c(auc_0_1$auc, auc_1_2$auc,auc_2_3$auc,auc_3_4$auc,auc_4_5$auc)), size=2) +
  geom_line(aes(x=0:4,y=apply(matrix_all_auc_bb, 2, mean), color="BB_all", linetype="all")) +
  geom_line(aes(x=0:4,y=apply(matrix_all_auc_bo, 2, mean), color="BO_all", linetype='all')) +
  geom_line(aes(x=0:4,y=apply(matrix_nodre_auc_bb, 2, mean), color="BB_nodre", linetype='nodre')) +
  geom_line(aes(x=0:4,y=apply(matrix_nodre_auc_bo, 2, mean), color="BO_nodre", linetype='nodre')) + 
  xlab('Follow up time (years)') + ylab("AUC") + ylim(0,1) +
  theme_bw()


ggplot() +
  geom_point(aes(x=0:4,y=c(prederr_0_1$prederr, prederr_1_2$prederr, prederr_2_3$prederr,
                           prederr_3_4$prederr, prederr_4_5$prederr)), size=2) +
  geom_line(aes(x=0:4,y=apply(matrix_all_prederr_bb, 2, mean), color="BB_all", linetype="all")) +
  geom_line(aes(x=0:4,y=apply(matrix_all_prederr_bo, 2, mean), color="BO_all", linetype='all')) +
  geom_line(aes(x=0:4,y=apply(matrix_nodre_prederr_bb, 2, mean), color="BB_nodre", linetype='nodre')) +
  geom_line(aes(x=0:4,y=apply(matrix_nodre_prederr_bo, 2, mean), color="BO_nodre", linetype='nodre')) + 
  xlab('Follow up time (years)') + ylab("AUC") +
  theme_bw()

