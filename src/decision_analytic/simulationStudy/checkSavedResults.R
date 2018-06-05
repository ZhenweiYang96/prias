load("Rdata/decision_analytic/PSA_Only/mvJoint_psa_light.Rdata")
savedFiles = list.files(path = "/home/a_tomer/Results/both_psa_dre_postMeans_mixed_censoring_tdist/", full.names = T)

postMeans = vector("list", length(savedFiles))
postwMeans = vector("list", length(savedFiles))
knots = vector("list", length(savedFiles))
coefficients = vector("list", length(savedFiles))

mcmcRes = vector("list", length(savedFiles))
weights = vector("list", length(savedFiles))
survModels = vector("list", length(savedFiles))

for(i in 1:length(savedFiles)){
  rm(jointModelData)
  load(savedFiles[i])
  print(paste("Reading the file number:", i))
  
  jointModelData$mvJoint_psa_simDs$mcmc$b=NULL
  #mcmcRes[[i]] = jointModelData$mvJoint_psa_simDs$mcmc
  #weights[[i]] = jointModelData$mvJoint_psa_simDs$weights
  
  postMeans[[i]] = jointModelData$mvJoint_psa_simDs$statistics$postMeans
  postwMeans[[i]] = jointModelData$mvJoint_psa_simDs$statistics$postwMeans
  knots[[i]] = jointModelData$mvJoint_psa_simDs$control$knots
  survModels[[i]] = jointModelData$survModel_simDs
  #coefficients[[i]] = jointModelData$jmFit$coefficients
  
  #postMeans[[i]] = jointModelData$mvJoint_dre_psa_simDs$statistics$postMeans
  #postwMeans[[i]] = jointModelData$mvJoint_dre_psa_simDs$statistics$postwMeans
  #knots[[i]] = jointModelData$mvJoint_dre_psa_simDs$control$knots
}

dev.off()
par(ask=T)
par(mfrow=c(2,2))
#Checking fixed effects
betas1 = sapply(postwMeans, function(x){x$betas1})
betas1 = sapply(coefficients, function(x){x$betas})
sapply(1:nrow(betas1), function(i){
  boxplot(betas1[i,])
  abline(h = mvJoint_psa_light$statistics$postwMeans$betas1[i])
})

sapply(1:nrow(betas1), function(i){
  temp = sapply(1:10, function(x){
    Hmisc::wtd.quantile(mcmcRes[[x]]$betas1[,i],weights = weights[[x]], type = "i/(n+1)", probs=0.5)
  })
  boxplot(temp)
  abline(h = mvJoint_psa_light$statistics$postwMeans$betas1[i])
})

#checking alphas
alphas = sapply(postwMeans, function(x){x$alphas})
alphas = rbind(sapply(coefficients, function(x){x$alpha}),sapply(coefficients, function(x){x$Dalpha}))
sapply(1:nrow(alphas), function(i){
  boxplot(alphas[i,])
  abline(h = mvJoint_psa_light$statistics$postwMeans$alphas[i])
})

sapply(1:nrow(alphas), function(i){
  temp = sapply(1:10, function(x){
    #Hmisc::wtd.quantile(mcmcRes[[x]]$alphas[,i],weights = weights[[x]], type = "i/(n+1)", probs=0.5)
    Hmisc::wtd.mean(mcmcRes[[x]]$alphas[,i],weights = weights[[x]])
  })
  boxplot(temp)
  abline(h = mvJoint_psa_light$statistics$postwMeans$alphas[i])
})


#checking gammas
gammas = sapply(postwMeans, function(x){x$gammas})
sapply(1:nrow(gammas), function(i){
  boxplot(gammas[i,])
  abline(h = mvJoint_psa_light$statistics$postMeans$gammas[i])
})

#checking variance covariance
sigma1 = sapply(postwMeans, function(x){x$sigma1})
boxplot(sigma1)
abline(h = mvJoint_psa_light$statistics$postMeans$sigma1)

D = sapply(postwMeans, function(x){c(x$D)})
sapply(1:nrow(D), function(i){
  boxplot(D[i,])
  abline(h = c(mvJoint_psa_light$statistics$postMeans$D)[i])
})
