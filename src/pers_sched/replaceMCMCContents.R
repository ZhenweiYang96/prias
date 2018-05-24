replaceMCMCContents <- function(fromObj, toObj){
  if(inherits(fromObj, "mvJMbayes") & inherits(toObj, "JMbayes")){
    #alpha
    alphas <- as.matrix(fromObj$mcmc[["alphas"]][,1])
    colnames(alphas) <- colnames(fromObj$mcmc[["alphas"]])[1]
    #colnames(alphas) <- "Assoct"
    
    if(ncol(fromObj$mcmc[["alphas"]]) > 1){
      Dalphas <- as.matrix(fromObj$mcmc[["alphas"]][,2:ncol(fromObj$mcmc[["alphas"]])])
      colnames(Dalphas) <- colnames(fromObj$mcmc[["alphas"]])[2:ncol(fromObj$mcmc[["alphas"]])]
      #colnames(Dalphas) <- "AssoctE"
    }
    
    #sigma
    sigma <- as.matrix(fromObj$mcmc[["sigma1"]])
    colnames(sigma) <- "sigma"
    
    #gammas
    gammas <- as.matrix(fromObj$mcmc[["gammas"]])
    
    #betas
    betas <- as.matrix(fromObj$mcmc[["betas1"]])
    
    #Bs.gammas
    Bs.gammas <- fromObj$mcmc[["Bs_gammas"]]
    colnames(Bs.gammas) <- paste("Bs.gammas", 1:ncol(Bs.gammas), sep="")
    
    #tauBs
    tauBs <- as.matrix(fromObj$mcmc[["tau_Bs_gammas"]])
    colnames(tauBs) <- "tauBs"
    
    #b matrix
    b <- fromObj$mcmc[["b"]]
    
    #D matrix
    dimD <- dim(fromObj$mcmc[["D"]])[3]
    D <- as.matrix(fromObj$mcmc[["D"]][,,1])
    if(dimD > 1){
      for(i in 2:(dimD)){
        D <- cbind(D, fromObj$mcmc[["D"]][,,i]) 
      }
    }
    
    for(i in 1:dimD){
      for(j in 1:dimD){
        colnames(D)[(i-1)*dimD + j] <- paste("D[", j, ", ", i, "]", sep="")
      }
    }
    
    #MCMC Obj
    toObj$mcmc[["betas"]] = betas
    toObj$mcmc[["sigma"]] = sigma 
    toObj$mcmc[["b"]] = b
    toObj$mcmc[["D"]] = D
    toObj$mcmc[["gammas"]] = gammas
    toObj$mcmc[["Bs.gammas"]] = Bs.gammas
    toObj$mcmc[["tauBs"]] = tauBs
    toObj$mcmc[["alphas"]] = alphas
    if(ncol(fromObj$mcmc[["alphas"]]) > 1){
      toObj$mcmc[["Dalphas"]] = Dalphas
    }

    toObj$postMeans <- lapply(toObj$mcmc, function (x) {
      d <- dim(x)
      if (!is.null(d) && length(d) > 2) apply(x, c(1L, 2L), mean) else colMeans(as.matrix(x))
    })
    dim(toObj$postMeans$D) <- c(dimD, dimD)
    colnames(toObj$postMeans$D) <- as.character(colnames(fromObj$statistics$postMeans$D))
    rownames(toObj$postMeans$D) <- colnames(toObj$postMeans$D)
    colnames(toObj$postMeans$b) <- colnames(toObj$postMeans$D)
    
    #Post modes
    toObj$postModes <- lapply(toObj$mcmc, function (x) apply(as.matrix(x), 2L, JMbayes:::modes))
    dim(toObj$postModes$D) <- c(dimD, dimD)
    
    #PostVarsRE
    toObj$postVarsRE <- apply(b, 1L, function (x) var(t(x)))
    dim(toObj$postVarsRE)<- c(dimD, dimD, dim(b)[[1]])
    
    #StErr
    #toObj$StErr = lapply(toObj$mcmc, JMbayes:::stdErr)
    
    #EffectiveSize
    #toObj$EffectiveSize = lapply(toObj$mcmc, JMbayes:::effectiveSize)
    
    #StDev
    #toObj$StDev = lapply(toObj$mcmc, function (x) apply(as.matrix(x), 2L, sd))
    
    #CIs
    #toObj$CIs = lapply(toObj$mcmc, function (x) 
     # apply(as.matrix(x), 2L, quantile, probs = c(0.025, 0.975)))
    
    #Pvalues
    #toObj$Pvalues = lapply(toObj$mcmc, function (x) apply(as.matrix(x), 2L, JMbayes:::computeP))
    
    return(toObj)
  }else{
    stop("The 'fromObj' parameter should be an mvJMbayes object and the 'toObj' parameter should be a JMbayes object.")
  }
}