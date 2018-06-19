load("Rdata/decision_analytic/PSA_Only/mvJoint_psa_easyrand_superlight.Rdata")
source("src/decision_analytic/load_lib.R")
source("src/decision_analytic/simulationStudy/timeDepSimDs_onlypsa.R")

savedFiles = list.files(path = "Rdata/decision_analytic/Simulation/Mixed/", full.names = T, pattern = "Rdata")
totalSavedFiles = length(savedFiles)

postMeans = vector("list", totalSavedFiles)
postwMeans = vector("list", totalSavedFiles)
jmFitCoefficients = vector("list", totalSavedFiles)
lmeFitCoefficients = vector("list", totalSavedFiles)
mvglmerFitCoefficients = vector("list", totalSavedFiles)

trueTimeDepCoeff = matrix(nrow = 4, ncol = totalSavedFiles)
postMeanTimeDepCoeff = matrix(nrow = 4, ncol = totalSavedFiles)
postwMeanTimeDepCoeff = matrix(nrow = 4, ncol = totalSavedFiles)
mvglmerTimeDepCoeff = matrix(nrow = 4, ncol = totalSavedFiles)
jmTimeDepCoeff = matrix(nrow=4, ncol=totalSavedFiles)
lmeTimeDepCoeff = matrix(nrow=4, ncol=totalSavedFiles)

for(i in 1:totalSavedFiles){
  rm(jointModelData)
  rm(timeDepDs)
  load(savedFiles[i])
  print(paste("Reading the file number:", i))
  
  postMeans[[i]] = jointModelData$mvJoint_psa_simDs$statistics$postMeans
  postMeans[[i]]$b = NULL
  
  postwMeans[[i]] = jointModelData$mvJoint_psa_simDs$statistics$postwMeans
  postwMeans[[i]]$b = NULL
  
  jmFitCoefficients[[i]] = jointModelData$jmFit$coefficients
  
  mvglmerFitCoefficients[[i]] = jointModelData$mvglmer_psa_simDs$postMeans
  mvglmerFitCoefficients[[i]]$b = NULL
  
  lmeFitCoefficients[[i]] = list(fixed=jointModelData$lmeFit$coefficients$fixed, 
                                 sigma=jointModelData$lmeFit$sigma,
                                 D = getVarCov(jointModelData$lmeFit))
  
  print("extracting time dependent ds results")
  #make time depDS
  timeDepDs = makeTimeDepDs(jointModelData$trainingData$trainingDs.id)
  trueTimeDepCoeff[, i] = getTimeDepDsCoefficients(timeDepDs, mvJoint_psa_easyrand_superlight$statistics$postwMeans$betas1)
  
  jointModelData$trainingData$trainingDs.id[,c("b_Int_PSA", "b_Slope1_PSA", "b_Slope2_PSA")] =  jointModelData$mvJoint_psa_simDs$statistics$postMeans$b
  rm(timeDepDs)
  timeDepDs = makeTimeDepDs(jointModelData$trainingData$trainingDs.id)
  postMeanTimeDepCoeff[, i] = getTimeDepDsCoefficients(timeDepDs, jointModelData$mvJoint_psa_simDs$statistics$postMeans$betas1)
  
  jointModelData$trainingData$trainingDs.id[,c("b_Int_PSA", "b_Slope1_PSA", "b_Slope2_PSA")] =  jointModelData$mvJoint_psa_simDs$statistics$postwMeans$b
  rm(timeDepDs)
  timeDepDs = makeTimeDepDs(jointModelData$trainingData$trainingDs.id)
  postwMeanTimeDepCoeff[, i] = getTimeDepDsCoefficients(timeDepDs, jointModelData$mvJoint_psa_simDs$statistics$postwMeans$betas1)

  jointModelData$trainingData$trainingDs.id[,c("b_Int_PSA", "b_Slope1_PSA", "b_Slope2_PSA")] =  jointModelData$mvglmer_psa_simDs$postMeans$b
  rm(timeDepDs)
  timeDepDs = makeTimeDepDs(jointModelData$trainingData$trainingDs.id)
  mvglmerTimeDepCoeff[,i] = getTimeDepDsCoefficients(timeDepDs, jointModelData$mvglmer_psa_simDs$postMeans$betas1)
    
  jointModelData$trainingData$trainingDs.id[,c("b_Int_PSA", "b_Slope1_PSA", "b_Slope2_PSA")] =  jointModelData$jmFit$EB$post.b
  rm(timeDepDs)
  timeDepDs = makeTimeDepDs(jointModelData$trainingData$trainingDs.id)
  jmTimeDepCoeff[, i] = getTimeDepDsCoefficients(timeDepDs, jointModelData$jmFit$coefficients$betas)
  
  jointModelData$trainingData$trainingDs.id[,c("b_Int_PSA", "b_Slope1_PSA", "b_Slope2_PSA")] =  random.effects(jointModelData$lmeFit)
  rm(timeDepDs)
  timeDepDs = makeTimeDepDs(jointModelData$trainingData$trainingDs.id)
  lmeTimeDepCoeff[, i] = getTimeDepDsCoefficients(timeDepDs, jointModelData$lmeFit$coefficients$fixed)
}

dev.off()
par(ask=T)

#Checking fixed effects
for(j in 1:length(postMeans[[1]]$betas1)){
  trueValue = mvJoint_psa_easyrand_superlight$statistics$postwMeans$betas1[j]
  
  postBeta = sapply(postMeans, function(x){x$betas1})[j, ]
  postwBeta = sapply(postwMeans, function(x){x$betas1})[j, ]
  postmvGlmerBeta = sapply(mvglmerFitCoefficients, function(x){x$betas1})[j, ]
  lmeBeta = sapply(lmeFitCoefficients, function(x){x$fixed})[j, ]
  jmFitBeta = sapply(jmFitCoefficients, function(x){x$betas})[j, ]
  
  plotDf = data.frame(beta = c(postBeta, postwBeta, postmvGlmerBeta, lmeBeta, jmFitBeta),
             Method = rep(c("postMean", "postwMean", "mvglmer_postMean", "lme", "JM"), each=totalSavedFiles))
  
  graph = ggplot(data=plotDf) + geom_boxplot(aes(x = Method, y=beta)) + ylab(attributes(postMeans[[1]]$betas1[j])$names[1]) +
    geom_hline(yintercept = trueValue, color="red")
  
  print(graph)
}

#checking alphas
for(j in 1:length(postMeans[[1]]$alphas)){
  trueValue = mvJoint_psa_easyrand_superlight$statistics$postwMeans$alphas[j]
  
  postAlpha = sapply(postMeans, function(x){x$alphas})[j, ]
  postwAlpha = sapply(postwMeans, function(x){x$alphas})[j, ]
  jmFitAlpha = sapply(jmFitCoefficients, function(x){c(x$alpha, x$Dalpha)})[j, ]
  
  trueTimeDepDsAlpha = trueTimeDepCoeff[j+2,]
  postMeanTimeDepDsAlpha = postMeanTimeDepCoeff[j+2,]
  postwMeanTimeDepDsAlpha = postwMeanTimeDepCoeff[j+2,]
  mvglmerTimeDepDsAlpha = mvglmerTimeDepCoeff[j+2,]
  jmTimeDepDsAlpha = jmTimeDepCoeff[j+2, ]
  lmeTimeDepDsAlpha = lmeTimeDepCoeff[j+2, ]
  
  plotDf = data.frame(alpha = c(postAlpha, postwAlpha, jmFitAlpha, 
                                trueTimeDepDsAlpha, postMeanTimeDepDsAlpha,mvglmerTimeDepDsAlpha,
                                jmTimeDepDsAlpha, lmeTimeDepDsAlpha),
                      color = rep(c("Bayes", "Bayes", "Freq", "timeDep", "timeDep",
                                    "timeDep", "timeDep", "timeDep"), each=totalSavedFiles),
                      Method = rep(c("postMean", "postwMean", "JM",
                                     "trueTimeDepDsAlpha", "postMeanTimeDepDsAlpha",
                                     "mvglmerTimeDepDsAlpha","jmTimeDepDsAlpha", "lmeTimeDepDsAlpha"),
                                   each=totalSavedFiles))
  
  graph = ggplot(data=plotDf) + geom_boxplot(aes(x = Method, y=alpha, fill=color)) + ylab(attributes(postMeans[[1]]$alphas[j])$names[1]) +
    geom_hline(yintercept = trueValue, color="red")
  
  print(graph)
}

#checking gammas
for(j in 1:length(postMeans[[1]]$gammas)){
  trueValue = mvJoint_psa_easyrand_superlight$statistics$postwMeans$gammas[j]
  
  postGamma = sapply(postMeans, function(x){x$gammas})[j, ]
  postwGamma = sapply(postwMeans, function(x){x$gammas})[j, ]
  jmFitGamma = sapply(jmFitCoefficients, function(x){x$gammas[2:3]})[j, ]
  
  trueTimeDepDsGamma = trueTimeDepCoeff[j+2,]
  postMeanTimeDepDsGamma = postMeanTimeDepCoeff[j+2,]
  postwMeanTimeDepDsGamma = postwMeanTimeDepCoeff[j+2,]
  mvglmerTimeDepDsGamma = mvglmerTimeDepCoeff[j+2,]
  jmTimeDepDsGamma = jmTimeDepCoeff[j+2, ]
  lmeTimeDepDsGamma = lmeTimeDepCoeff[j+2, ]
  
  plotDf = data.frame(gamma = c(postGamma, postwGamma, jmFitGamma, 
                                trueTimeDepDsGamma, postMeanTimeDepDsGamma,mvglmerTimeDepDsGamma,
                                jmTimeDepDsGamma, lmeTimeDepDsGamma),
                      color = rep(c("Bayes", "Bayes", "Freq", "timeDep", "timeDep",
                                    "timeDep", "timeDep", "timeDep"), each=totalSavedFiles),
                      Method = rep(c("postMean", "postwMean", "JM",
                                     "trueTimeDepDsGamma", "postMeanTimeDepDsGamma",
                                     "mvglmerTimeDepDsGamma","jmTimeDepDsGamma", "lmeTimeDepDsGamma"),
                                   each=totalSavedFiles))
  
  graph = ggplot(data=plotDf) + geom_boxplot(aes(x = Method, y=gamma, fill=color)) + ylab(attributes(postMeans[[1]]$gammas[j])$names[1]) +
    geom_hline(yintercept = trueValue, color="red")
  
  print(graph)
}

#checking sigma
trueValue = mvJoint_psa_easyrand_superlight$statistics$postwMeans$sigma1

postSigma= sapply(postMeans, function(x){x$sigma1})
postwSigma = sapply(postwMeans, function(x){x$sigma1})
postmvGlmerSigma = sapply(mvglmerFitCoefficients, function(x){x$sigma1})
lmeSigma = sapply(lmeFitCoefficients, function(x){x$sigma})
jmFitSigma = sapply(jmFitCoefficients, function(x){x$sigma})

plotDf = data.frame(sigma = c(postSigma, postwSigma, postmvGlmerSigma, lmeSigma, jmFitSigma),
                    Method = rep(c("postMean", "postwMean", "mvglmer_postMean", "lme", "JM"), each=totalSavedFiles))

graph = ggplot(data=plotDf) + geom_boxplot(aes(x = Method, y=sigma)) + ylab("Sigma") +
  geom_hline(yintercept = trueValue, color="red")
print(graph)

#Checking D matrix
for(j in 1:nrow(postMeans[[1]]$D)){
  for(k in j:ncol(postMeans[[1]]$D)){
    trueValue = mvJoint_psa_easyrand_superlight$statistics$postwMeans$D[j,k]
    
    postD = sapply(postMeans, function(x){x$D[j,k]})
    postwD = sapply(postwMeans, function(x){x$D[j,k]})
    postmvGlmerD = sapply(mvglmerFitCoefficients, function(x){x$D[j,k]})
    lmeD = sapply(lmeFitCoefficients, function(x){x$D[j,k]})
    jmFitD = sapply(jmFitCoefficients, function(x){x$D[j,k]})
    
    plotDf = data.frame(D = c(postD, postwD, postmvGlmerD, lmeD, jmFitD),
                        Method = rep(c("postMean", "postwMean", "mvglmer_postMean", "lme", "JM"), each=totalSavedFiles))
    
    graph = ggplot(data=plotDf) + geom_boxplot(aes(x = Method, y=D)) + ylab(paste0("D_", j,"_", k)) +
      geom_hline(yintercept = trueValue, color="red")
    
    print(graph)
  }
}