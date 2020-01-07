files = list.files("Rdata/lastpaper/auc_pe/",
                   pattern = "auc", full.names = T)

t_horizs = seq(1, 6, 0.5)

aucs = lapply(files, FUN = function(file){
  load(file)
  return(sapply(auc_list, "[[", "auc"))
})

aucs = do.call('rbind', aucs)
colnames(aucs) = t_horizs
apply(aucs, 2, FUN = mean)
apply(aucs, 2, FUN = quantile, probs=0.025)
apply(aucs, 2, FUN = quantile, probs=0.975)

write.csv(round(t(rbind(apply(aucs, 2, FUN = mean),
                        apply(aucs, 2, FUN = quantile, probs=0.025),
                        apply(aucs, 2, FUN = quantile, probs=0.975))) ,3),
          row.names = T, file = "temp.csv")

######################
files = list.files("Rdata/lastpaper/auc_pe/",
                   pattern = "pe", full.names = T)

t_horizs = seq(1, 6, 0.5)

mape = lapply(files, FUN = function(file){
  load(file)
  return(sapply(sapply(pe_list, "[[", "ape"), mean, na.rm=T))
})

mape = do.call('rbind', mape)
colnames(mape) = t_horizs
apply(mape, 2, FUN = mean)
apply(mape, 2, FUN = quantile, probs=0.025)
apply(mape, 2, FUN = quantile, probs=0.975)

write.csv(round(t(rbind(apply(mape, 2, FUN = mean),
                        apply(mape, 2, FUN = quantile, probs=0.025),
                        apply(mape, 2, FUN = quantile, probs=0.975))) ,3),
          row.names = T, file = "temp.csv")
