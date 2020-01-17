dim(auc[[2]]$UCSF)

tt=data.frame("Recalibrated", "UCSF",
              as.numeric(colnames(auc[[2]]$UCSF))-1,
              as.numeric(colnames(auc[[2]]$UCSF)),
              apply(auc[[2]]$UCSF, MARGIN = 2, mean, na.rm=T),
              apply(auc[[2]]$UCSF, MARGIN = 2, quantile, probs=0.025, na.rm=T),
              apply(auc[[2]]$UCSF, MARGIN = 2, quantile, probs=0.975, na.rm=T))
colnames(tt) = colnames(auc_df)

auc_df2 = rbind(auc_df, tt)