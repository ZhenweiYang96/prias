temp = dataset.id$`1`
temp = droplevels(temp)

survfit_func = function(id, times){
  survfitJM(fitted_JM, dataset[[as.character(id)]], last.time = 0, survTimes = times, idVar="P_ID")[[1]][[1]][,"Mean"]
}

temp$auc = NA
for(i in 1:nrow(temp)){
  temp$auc[i] = integrate(survfit_func, lower = 1, upper = 5, id =temp$P_ID[i])$value
}

temp$surv1 = temp$surv2= NA
for(i in 1:nrow(temp)){
  temp[i, c('surv1', 'surv2')] = survfitJM(fitted_JM,dataset[[as.character(temp$P_ID[i])]], last.time = 0, survTimes = c(1,2), idVar="P_ID")[[1]][[1]][,"Mean"]
}
