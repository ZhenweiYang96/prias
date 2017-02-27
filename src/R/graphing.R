idList = droplevels(unique(amctx.id[amctx.id$gl_failure==0,]$amctx))

idList = c(73, 94, 195, 209)
for(id in idList){
  ND = amctx_creatinine[amctx_creatinine$amctx==id,]
  futureTimes = seq(max(ND$tx_s_years), (max(ND$tx_s_years) + 3), 0.1)
  
  #sfit.patient2 = survfitJM(jointfit_creatinine_tdboth_nomv, ND, idVar="amctx", survTimes = futureTimes)
  #plot(sfit.patient2, estimator="mean", include.y=T, conf.int=T, fill.area=T, col.area="lightgrey", main=paste("amctx =",id))
  #Sys.sleep(2)
  
  longprof = predict(jointfit_creatinine_tdboth_nomv, ND, type = "Subject",
          interval = "confidence", return = TRUE, idVar="amctx", FtTimes = futureTimes)
  last.time <- with(longprof, tx_s_years[!is.na(low)][1])
  lattice::xyplot(pred + low + upp ~ tx_s_years, data = longprof, type = "l", 
                  lty = c(1, 2, 2), col = c(2, 1, 1), abline = list(v = last.time, lty = 3),
            xlab = "Time (years)", ylab = "Predicted log(serum creatinine)", main=paste("amctx =",id))
}

