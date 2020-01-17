# ggplot(mape_df) + geom_line(aes(x=t_horiz, y=mean, color=cohort)) +
#   geom_ribbon(aes(x=t_horiz, ymin=lower, ymax=upper), alpha=0.5) +
#   facet_grid(model~cohort) + ylim(0,1) + xlim(0,10)

setwd("C:/Users/838035/Google Drive/PhD/src/prias")
pe_bootstrap_df = do.call('rbind', lapply(cohortnames, function(cohort){
  t_horizs = seq(1, round(reclassification_df$time_10pat_risk_set[reclassification_df$Cohort==cohort]), 0.5)

  cohort_df = lapply(1:30, function(iter){
    seed = 2019 + iter
    load(paste0("Rdata/gap3/PRIAS_2019/validation/pe/pe_gof_recalib_model/", cohort, "_", seed, ".Rdata"))
    for(i in 1:length(t_horizs)){
      pe_list[[i]] = pe_list[[i]][, c("P_ID", "right_cens_time", "reclassification", 
                                      "real_period_status","cum_risk_T_start_T_horiz", 
                                      "ape", "spe")]
      pe_list[[i]]$t_start = t_horizs[i] - 1
      pe_list[[i]]$t_horiz = t_horizs[i]
      pe_list[[i]]$bs_iter = iter
      pe_list[[i]]$center = cohort
    }
    return(do.call('rbind', pe_list))
  })

  return(do.call('rbind', cohort_df))
}))

pe_bootstrap_df = do.call('rbind', by(pe_bootstrap_df$center, data = pe_bootstrap_df, FUN = function(cohort_df){
  cohort_df_list = by(cohort_df$t_horiz, data=cohort_df, FUN = function(t_horiz_df){
    mape_t_horiz = by(t_horiz_df$bs_iter, data=t_horiz_df, FUN = function(bs_iter_df){
      mape = sum(bs_iter_df$ape)/nrow(bs_iter_df)
      mspe = sum(bs_iter_df$spe)/nrow(bs_iter_df)
      return(data.frame(cohort = bs_iter_df$center[1],
                        t_start = bs_iter_df$t_start[1],
                        t_horiz = bs_iter_df$t_horiz[1],
                        bs_iter = bs_iter_df$bs_iter[1],
                        mape = mape,
                        mspe = mspe))
    })
    mape_t_horiz = do.call('rbind', mape_t_horiz)
    return(mape_t_horiz)
  })
  return(do.call('rbind', cohort_df_list))
}))

save(pe_bootstrap_df, file="Rdata/gap3/PRIAS_2019/validation/pe_bootstrap_df.Rdata")

pe_df = do.call('rbind',by(pe_bootstrap_df$cohort, data=pe_bootstrap_df, FUN = function(x){
  return(do.call('rbind',by(x$t_horiz, data = x, FUN = function(y){
    z = y[1,1:3]
    z$mean_mape = mean(y$mape, na.rm = T)
    z$lower_mape = quantile(y$mape, probs = 0.025, na.rm = T)
    z$upper_mape = quantile(y$mape, probs = 0.975, na.rm = T)
    
    z$mean_mspe = mean(y$mspe,na.rm = T)
    z$lower_mspe = quantile(y$mspe, probs = 0.025,na.rm = T)
    z$upper_mspe = quantile(y$mspe, probs = 0.975,na.rm = T)
    return(z)
  })))
}))

save(pe_df, file="Rdata/gap3/PRIAS_2019/validation/pe_df.Rdata")
