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

#######################################################################
# the following function creates the predicted values
# and the 95% CIs
effectPlotData <- function (object, newdata, orig_data) {
  form <- formula(object)
  namesVars <- all.vars(form)
  betas <- if (!inherits(object, "lme")) coef(object) else fixef(object)
  V <- if (inherits(object, "geeglm")) object$geese$vbeta else vcov(object)
  orig_data <- orig_data[complete.cases(orig_data[namesVars]), ]
  Terms <- delete.response(terms(form))
  mfX <- model.frame(Terms, data = orig_data)
  Terms_new <- attr(mfX, "terms")
  mfX_new <- model.frame(Terms_new, newdata, xlev = .getXlevels(Terms, mfX))
  X <- model.matrix(Terms_new, mfX_new)
  pred <- c(X %*% betas)
  ses <- sqrt(diag(X %*% V %*% t(X)))
  newdata$pred <- pred
  newdata$low <- pred - 1.96 * ses
  newdata$upp <- pred + 1.96 * ses
  newdata
}

# the data frame that contains the combination of values to
# create the plot
newDF <- with(amctx, expand.grid(tx_s_years = seq(0, 15, length.out = 30),
                                rec_gender = "M", 
                                rec_age_fwp1 = median(amctx.id$rec_age), 
                                d_age = median(amctx.id$d_age),
                                d_bmi = median(amctx.id$d_bmi),
                                rec_bmi = median(amctx.id$rec_bmi),
                                tx_dial_days = median(amctx.id$tx_dial_days),
                                tx_cit = median(amctx.id$tx_cit),
                                tx_pra = median(amctx.id$tx_pra),
                                tx_dm = "no", tx_previoustx = "no",
                                tx_hla = "6", is_cni = "yes", tx_dgf = "no"))

# the effects plot
xyplot(pred + low + upp ~ tx_s_years | rec_gender, 
       data = effectPlotData(model_creatinine, newDF, amctx), 
       lty = c(1, 2, 2), col = c(2, 1, 1), lwd = 2, type = "l",
       xlab = "Follow-up time (years)",
       ylab = "log (serum creatinine)")