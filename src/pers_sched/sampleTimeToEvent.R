rConditionalTimeToEvent = function (object, newdata, idVar = "id", simulate = TRUE,
                              last.time = NULL, LeftTrunc_var = NULL, M = 200L, 
                              scale = 1.6, init.b = NULL, seed = 1L, ...) {
  if (!inherits(object, "JMbayes"))
    stop("Use only with 'JMbayes' objects.\n")
  if (!is.data.frame(newdata) || nrow(newdata) == 0L)
    stop("'newdata' must be a data.frame with more than one rows.\n")
  if (is.null(newdata[[idVar]]))
    stop("'idVar' not in 'newdata.\n'")
  
  TT <- object$y$Time
  timeVar <- object$timeVar
  df.RE <- object$y$df.RE
  param <- object$param
  densLong <- object$Funs$densLong
  hasScale <- object$Funs$hasScale
  anyLeftTrunc <- object$y$anyLeftTrunc
  densRE <- object$Funs$densRE
  transFun.value <- object$Funs$transFun.value
  transFun.extra <- object$Funs$transFun.extra
  extraForm <- object$Forms$extraForm
  indFixed <- extraForm$indFixed
  indRandom <- extraForm$indRandom
  indBetas <- object$y$indBetas
  TermsX <- object$Terms$termsYx
  TermsZ <- object$Terms$termsYz
  TermsX.extra <- object$Terms$termsYx.extra
  TermsZ.extra <- object$Terms$termsYz.extra
  mfX <- model.frame.default(TermsX, data = newdata)
  mfZ <- model.frame.default(TermsZ, data = newdata)
  formYx <- reformulate(attr(delete.response(TermsX), "term.labels"))
  formYz <- object$Forms$formYz
  estimateWeightFun <- object$estimateWeightFun
  weightFun <- object$Funs$weightFun
  max.time <- max(TT)
  na.ind <- as.vector(attr(mfX, "na.action"))
  na.ind <- if (is.null(na.ind)) {
    rep(TRUE, nrow(newdata))
  } else {
    !seq_len(nrow(newdata)) %in% na.ind
  }
  id <- as.numeric(unclass(newdata[[idVar]]))
  id <- id. <- match(id, unique(id))
  id <- id[na.ind]
  y <- model.response(mfX)
  X <- model.matrix.default(formYx, mfX)
  Z <- model.matrix.default(formYz, mfZ)[na.ind, , drop = FALSE]
  TermsT <- object$Terms$termsT
  data.id <- newdata[tapply(row.names(newdata), id, tail, n = 1L), ]
  data.s <- data.id[rep(1:nrow(data.id), each = object$control$GQsurv.k), ]
  idT <- data.id[[idVar]]
  idT <- match(idT, unique(idT))
  ids <- data.s[[idVar]]
  ids <- match(ids, unique(ids))
  
  mfT <- model.frame.default(delete.response(TermsT), data = data.id)
  formT <- if (!is.null(kk <- attr(TermsT, "specials")$cluster)) {
    tt <- drop.terms(TermsT, kk - 1, keep.response = FALSE)
    reformulate(attr(tt, "term.labels"))
  } else {
    tt <- attr(delete.response(TermsT), "term.labels")
    if (length(tt)) reformulate(tt) else reformulate("1")
  }
  
  W <- model.matrix.default(formT, mfT)[, -1L, drop = FALSE]
  obs.times <- split(newdata[[timeVar]][na.ind], id)
  last.time <- if (is.null(last.time)) {
    tapply(newdata[[timeVar]], id., tail, n = 1L)
  } else if (is.character(last.time) && length(last.time) == 1L) {
    tapply(newdata[[last.time]], id., tail, n = 1L)
  } else if (is.numeric(last.time)) {
    rep_len(last.time, length.out = nrow(data.id))
  } else {
    stop("\nnot appropriate value for 'last.time' argument.")
  }
  
  TimeL <- if (!is.null(anyLeftTrunc) && anyLeftTrunc) {
    if (is.null(LeftTrunc_var) || is.null(newdata[[LeftTrunc_var]])) {
      warning("The original joint model was fitted in a data set with left-",
              "truncation and\nargument 'LeftTrunc_var' of survfitJM() has not ", 
              "been specified.\n")
    }
    TimeL <- newdata[[LeftTrunc_var]]
    tapply(TimeL, id, head, n = 1)
  }
  
  n <- length(TT)
  n.tp <- length(last.time)
  ncx <- ncol(X)
  ncz <- ncol(Z)
  ncww <- ncol(W)
  if (ncww == 0L)
    W <- NULL
  lag <- object$y$lag
  betas <- object$postMeans$betas
  sigma <- object$postMeans$sigma
  D <- object$postMeans$D
  gammas <- object$postMeans$gammas
  alphas <- object$postMeans$alphas
  Dalphas <- object$postMeans$Dalphas
  shapes <- object$postMeans$shapes
  Bs.gammas <- object$postMeans$Bs.gammas
  list.thetas <- list(betas = betas, sigma = sigma, gammas = gammas, alphas = alphas, 
                      Dalphas = Dalphas, shapes = shapes, Bs.gammas = Bs.gammas, D = D)
  list.thetas <- list.thetas[!sapply(list.thetas, is.null)]
  thetas <- unlist(as.relistable(list.thetas))
  
  #######################################
  # To be removed later
  #######################################
  logh.b = JMbayes:::logh.b
  environment(logh.b) <- asNamespace('JMbayes')
  S.b = JMbayes:::S.b
  environment(S.b) <- asNamespace('JMbayes')
  log.posterior.b = JMbayes:::log.posterior.b
  environment(log.posterior.b) <- asNamespace('JMbayes')
  hMats = JMbayes:::hMats
  environment(hMats) <- asNamespace('JMbayes')
  ModelMats = JMbayes:::ModelMats
  environment(ModelMats) <- asNamespace('JMbayes')
  gaussKronrod = JMbayes:::gaussKronrod
  environment(gaussKronrod) <- asNamespace('JMbayes')
  fastSumID = JMbayes:::fastSumID
  environment(fastSumID) <- asNamespace('JMbayes')
  rmvt = JMbayes:::rmvt
  environment(rmvt) <- asNamespace('JMbayes')
  dmvt = JMbayes:::dmvt
  environment(dmvt) <- asNamespace('JMbayes')
  ###############################################

  environment(log.posterior.b) <- environment(S.b) <- environment(logh.b) <- environment()
  environment(hMats) <- environment(ModelMats) <- environment()
  # construct model matrices to calculate the survival functions
  obs.times.surv <- split(data.id[[timeVar]], idT)
  survMats.last <- vector("list", n.tp)
  for (i in seq_len(n.tp)) {
    survMats.last[[i]] <- ModelMats(last.time[i], ii = i, timeL = TimeL[i])
  }
  
  # calculate the Empirical Bayes estimates and their (scaled) variance
  modes.b <- matrix(0, n.tp, ncz)
  invVars.b <- Vars.b <- vector("list", n.tp)
  for (i in seq_len(n.tp)) {
    betas.new <- betas
    sigma.new <- sigma
    D.new <- D
    gammas.new <- gammas
    alphas.new <- alphas
    Dalphas.new <- Dalphas
    shapes.new <- shapes
    Bs.gammas.new <- Bs.gammas
    ff <- function (b, y, tt, mm, i) -log.posterior.b(b, y, Mats = tt, ii = i)
    start <- if (is.null(init.b)) rep(0, ncz) else init.b[i, ]
    opt <- try(optim(start, ff, y = y, tt = survMats.last, i = i, 
                     method = "BFGS", hessian = TRUE), silent = TRUE)
    if (inherits(opt, "try-error")) {
      gg <- function (b, y, tt, mm, i) cd(b, ff, y = y, tt = tt, i = i)
      opt <- optim(start, ff, gg, y = y, tt = survMats.last, 
                   i = i, method = "BFGS", hessian = TRUE, 
                   control = list(parscale = rep(0.1, ncz)))
    } 
    modes.b[i, ] <- opt$par
    invVars.b[[i]] <- opt$hessian/scale
    Vars.b[[i]] <- scale * solve(opt$hessian)
  }
  
  set.seed(seed)
  out <- vector("list", M)
  success.rate <- matrix(FALSE, M, n.tp)
  b.old <- b.new <- modes.b
  old_Tj <- Tj <- 1.1 * last.time
  if (n.tp == 1)
    dim(b.old) <- dim(b.new) <- c(1L, ncz)
  mcmc <- object$mcmc
  mcmc <- mcmc[names(mcmc) != "b"]
  if (M > nrow(mcmc$betas)) {
    warning("'M' cannot be set greater than ", nrow(mcmc$betas))
    M <- nrow(mcmc$betas)
    out <- vector("list", M)
    success.rate <- matrix(FALSE, M, n.tp)
  }
  samples <- sample(nrow(mcmc$betas), M)
  mcmc[] <- lapply(mcmc, function (x) x[samples, , drop = FALSE])
  proposed.b <- mapply(rmvt, mu = split(modes.b, row(modes.b)), Sigma = Vars.b, 
                       MoreArgs = list(n = M, df = 4), SIMPLIFY = FALSE)
  proposed.b[] <- lapply(proposed.b, function (x) if (is.matrix(x)) x else rbind(x))
  dmvt.proposed <- mapply(dmvt, x = proposed.b, mu = split(modes.b, row(modes.b)),
                          Sigma = Vars.b, MoreArgs = list(df = 4, log = TRUE), 
                          SIMPLIFY = FALSE)
  
  log.p_Tj <- function (patientNum, Tj) {
    log.S_ti <- S.b(last.time[patientNum], b.new[patientNum, ] , i = patientNum, survMats.last[[patientNum]], log = TRUE)
    log.S_Tj <- S.b(Tj, b.new[patientNum, ] , i = patientNum, ModelMats(Tj, patientNum), log = TRUE)
    log.h_Tj <- logh.b(b.new[patientNum, ] , hMats(Tj))
    log.h_Tj + log.S_Tj - log.S_ti
  }
  
  sliceTj <- function (patientNum, current.Tj, step=1, maxR = 20){
    logy <- log.p_Tj(patientNum, current.Tj) - rexp(1L)
    tt <- runif(1L, 0, step)
    L <- current.Tj - tt
    R <- current.Tj + step - tt
    while (L > last.time[patientNum] && log.p_Tj(patientNum, L) > logy) {
      L <- L - step
    }
    while (log.p_Tj(patientNum, R) > logy & R<maxR) {
      R <- R + step
    }
    L <- max(last.time[patientNum], L)
    
    repeat {
      newTj <- runif(1L, L, R)
      new.logPost <- log.p_Tj(patientNum, newTj)
      if (new.logPost >= logy)
        break
      if (newTj < current.Tj) {
        L <- newTj
      }
      else {
        R <- newTj
      }
    }
    newTj
  }
  
  for (m in 1:M) {
    # Step 1: extract parameter values
    betas.new <- mcmc$betas[m, ]
    if (hasScale)
      sigma.new <- mcmc$sigma[m, ]
    if (!is.null(W))
      gammas.new <- mcmc$gammas[m, ]
    if (param %in% c("td-value", "td-both", "shared-betasRE", "shared-RE")) 
      alphas.new <- mcmc$alpha[m, ]
    if (param %in% c("td-extra", "td-both"))
      Dalphas.new <- mcmc$Dalphas[m, ]
    if (estimateWeightFun)
      shapes.new <- mcmc$shapes[m, ]
    D.new <- mcmc$D[m, ]; dim(D.new) <- dim(D)
    Bs.gammas.new <- mcmc$Bs.gammas[m, ]
    
    SS <- rep(NA, n.tp)
    for (i in seq_len(n.tp)) {
      # Step 2: simulate new random effects values
      p.b <- proposed.b[[i]][m, ]
      dmvt.old <- dmvt(b.old[i, ], modes.b[i, ], invSigma = invVars.b[[i]], 
                       df = 4, log = TRUE)
      dmvt.prop <- dmvt.proposed[[i]][m]
      a <- min(exp(log.posterior.b(p.b, y, survMats.last, ii = i) + dmvt.old - 
                     log.posterior.b(b.old[i, ], y, survMats.last, ii = i) - dmvt.prop), 1)
      ind <- runif(1) <= a
      success.rate[m, i] <- ind
      if (!is.na(ind) && ind)
        b.new[i, ] <- p.b
      
      # # Step 3: Simulate T_j^* from [T_j^* | T_j > last.time, Y_j(s)]
      prop_Tj <- runif(1, last.time[i], max.time * 1.1)
      aa <- min(exp(log.p_Tj(i, prop_Tj) - log.p_Tj(i, old_Tj[i])), 1)
      ind <- runif(1) <= aa
      if (!is.na(ind) && ind) {
        Tj[i] <- prop_Tj
      }
      
      #Tj[i] = sliceTj(i, runif(1, last.time[i], 20), step = 3)
      #Tj[i] = sliceTj(i, Tj[i], step = 3)
      
      SS[i] = Tj[i]
    }
    old_Tj <- Tj
    b.old <- b.new
    out[[m]] <- SS
  }

  res <- vector("list", n.tp)
  for (i in seq_len(n.tp)) {
    res[[i]] = unlist(lapply(out, function(x){x[[i]]}))
  }
  
  names(res) <- unique(unclass(newdata[[idVar]]))
  class(res) <- "rConditionalTimeToEvent.JMbayes"
  res
}

environment(rConditionalTimeToEvent) <- asNamespace('JMbayes')

tStart = Sys.time()
#eTime = expectedCondFailureTime(ND)
#varTime = varianceCondFailureTime(ND)
tSample = rConditionalTimeToEvent(joint_psa_replaced, ND, idVar = "P_ID", M = 1000)
tEnd = Sys.time()

tEnd-tStart
