# pid = sample(prias.id.rightCens[prias.id.rightCens$progressed >0,]$P_ID, size = 1)
# last.time = getLastBiopsyTime(pid, ifelse(prias.id.rightCens$progressed[prias.id.rightCens$P_ID==pid]>0,2,1))
# x = survfitJM(joint_psa_replaced_prias, psa_data_set[psa_data_set$P_ID==pid,], idVar = "P_ID", last.time = last.time)
# plot(x, estimator = "mean", include.y = TRUE,conf.int = TRUE, fill.area = TRUE, col.area = "lightgrey")
# 
# ds = psa_data_set[psa_data_set$P_ID %in% prias.id.rightCens$P_ID[prias.id.rightCens$progressed > 0],]
# ds$P_ID = droplevels(ds$P_ID)
# 
# ord = order(c(by(ds$psa, ds$P_ID, function(x){max(x)})),decreasing = T)
# unique(ds$P_ID)[ord]


pid = sample(prias.id.rightCens[prias.id.rightCens$progressed >0,]$P_ID, size = 1)
last.time = getLastBiopsyTime(pid, ifelse(prias.id.rightCens$progressed[prias.id.rightCens$P_ID==pid]>0,2,1))
x = survfitJM(joint_psa_replaced_prias, psa_data_set[psa_data_set$P_ID==pid,], idVar = "P_ID", last.time = last.time)
plot(x, estimator = "mean", include.y = TRUE,conf.int = TRUE, fill.area = TRUE, col.area = "lightgrey")

ds = psa_data_set[psa_data_set$P_ID %in% prias.id.rightCens$P_ID[prias.id.rightCens$progressed > 0],]
ds$P_ID = droplevels(ds$P_ID)

ord = order(c(by(ds$psa, ds$P_ID, function(x){max(x)})),decreasing = T)
unique(ds$P_ID)[ord]
