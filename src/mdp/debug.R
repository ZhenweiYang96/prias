survfit_temp <- function (object, newdata, survTimes = NULL, idVar = "id", 
                          last.time = NULL, M = 200L, scale = 1.6, log = FALSE, 
                          CI.levels = c(0.025, 0.975), seed = 2019, ...) {
  if (!inherits(object, "mvJMbayes"))
    stop("Use only with 'mvJMbayes' objects.\n")
  if (!is.data.frame(newdata) || nrow(newdata) == 0L)
    stop("'newdata' must be a data.frame with more than one rows.\n")
  if (is.null(newdata[[idVar]]))
    stop("'idVar' not in 'newdata'.\n")
  timeVar <- object$model_info$timeVar
  TermsU <- object$model_info$coxph_components$TermsU
  control <- object$control
  families <- object$model_info$families
  fams <- sapply(families, "[[", "family")
  links <- sapply(families, "[[", "link")
  n_outcomes <- length(object$model_info$families)
  seq_n_outcomes <- seq_len(n_outcomes)
  componentsL <- object$model_info$mvglmer_components
  componentsS <- object$model_info$coxph_components
  TermsL <- componentsL[grep("Terms", names(componentsL), fixed = TRUE)]
  TermsL <- lapply(TermsL, function (x) {environment(x) <- parent.frame(); x})
  TermsFormulas_fixed <- object$model_info$coxph_components$TermsFormulas_fixed
  TermsFormulas_random <- object$model_info$coxph_components$TermsFormulas_random
  build_model_matrix <- function (Terms, data) {
    out <- vector("list", length(Terms))
    for (i in seq_along(Terms)) {
      MF <- model.frame(Terms[[i]], data = data, na.action = NULL)
      out[[i]] <- model.matrix(Terms[[i]], MF)
    }
    out
  }
  #######################################################################
  mfLna <- lapply(TermsL, FUN = model.frame.default, data = newdata)
  mfLnaNULL <- lapply(TermsL, FUN = model.frame.default, data = newdata, na.action = NULL)
  na.inds <- lapply(mfLna, FUN = function (x) as.vector(attr(x, "na.action")))
  na.inds <- lapply(na.inds, FUN = function (x) {
    if (is.null(x)) {
      rep(TRUE, nrow(newdata))
    } else {
      !seq_len(nrow(newdata)) %in% x
    }
  })
  na.inds <- na.inds[grep("TermsX", names(na.inds))]
  #######################################################################
  TermsS <- componentsS[grep("Terms", names(componentsS), fixed = TRUE)][[1]]
  formulasL <- lapply(TermsL, FUN = function (x) formula(x))
  formulasL2 <- lapply(TermsL, FUN = function (x) formula(delete.response(x))) 
  formulasS <- formula(delete.response(TermsS))
  allvarsL <- unique(unlist(lapply(formulasL, FUN = all.vars)))
  allvarsS <- all.vars(formulasS)
  allvarsLS <- unique(c(allvarsL, allvarsS))
  respVars <- unlist(componentsL[grep("respVar", names(componentsL), fixed = TRUE)], use.names = FALSE)
  id <- as.numeric(unclass(newdata[[idVar]]))
  id <- id. <- match(id, unique(id))
  obs.times <- split(newdata[[timeVar]], id)
  obs.times[] <- mapply(function (x, nams) {names(x) <- nams; x}, obs.times, 
                        split(row.names(newdata), id), SIMPLIFY = FALSE)
  mfL <- lapply(TermsL, FUN = model.frame.default, data = newdata, na.action = NULL)
  y <- lapply(mfL[grep("TermsX", names(mfL))], model.response)
  if (any(!allvarsLS %in% names(newdata)))
    stop("The following variable or variables: ", 
         paste(c(allvarsL, allvarsS)[which(!c(allvarsL, allvarsS) %in% names(newdata))], "\n", sep = " "),
         "should be included in 'newdata'.\n")
  idNewL <- factor(newdata[[idVar]], levels = unique(newdata[[idVar]]))
  n.NewL <- length(unique(idNewL))
  na.inds2 <- vector("list", n.NewL)
  for (i in 1:n.NewL) {
    for (j in 1:n_outcomes) {
      tmp <- split(na.inds[[j]], id)
      na.inds2[[i]][[j]] <- unname(unlist(tmp[[i]])) 
    }
  }
  last.time <- if (is.null(last.time)) {
    tapply(newdata[[timeVar]], idNewL, FUN = max, simplify = FALSE)
  } else {
    rep_len(last.time, length.out = length(unique(newdata[[idVar]])))
  }
  Time <- componentsS$Time
  if (is.null(survTimes) || !is.numeric(survTimes)) {
    survTimes <- seq(min(Time), max(Time) + 0.01, length.out = 35L)
  }
  times.to.pred_upper <- lapply(last.time, FUN = function (t) survTimes[survTimes > t])
  times.to.pred_lower <- mapply(FUN = function (t1, t2) as.numeric(c(t1, t2[-length(t2)])), 
                                last.time, times.to.pred_upper, SIMPLIFY = FALSE)
  GQsurv <- if (control$GQsurv == "GaussKronrod") gaussKronrod() else gaussLegendre(control$GQsurv.k)
  wk <- GQsurv$wk
  sk <- GQsurv$sk
  K <- length(sk)
  P <- lapply(last.time, FUN = function (x) x / 2)
  P1 <- mapply(FUN = function (x, y) (x + y) / 2, times.to.pred_upper, times.to.pred_lower, SIMPLIFY = FALSE)
  P2 <- mapply(FUN = function (x, y) (x - y) / 2, times.to.pred_upper, times.to.pred_lower, SIMPLIFY = FALSE)
  GK_points_postRE <- matrix(unlist(lapply(P, FUN = function (x) outer(x, sk + 1))), 
                             ncol = K, byrow = TRUE)
  GK_points_CumHaz <- mapply(FUN = function (x, y, sk) outer(y, sk) + x, P1, P2, 
                             SIMPLIFY = F, MoreArgs = list(sk = sk))
  idGK <- rep(seq_len(n.NewL), each = K)
  
  newdata[[idVar]] <- match(newdata[[idVar]], unique(newdata[[idVar]]))
  
  newdata.GK.postRE <- right_rows(newdata, newdata[[timeVar]], newdata[[idVar]], GK_points_postRE)
  newdata.GK.postRE[[timeVar]] <- c(t(GK_points_postRE))
  newdata.GK.postRE[[idVar]] <- match(newdata.GK.postRE[[idVar]], unique(newdata.GK.postRE[[idVar]]))
  
  newdata.GK.CumHaz <- vector("list", n.NewL)
  for(i in 1:n.NewL) {
    upred.lngth <- nrow(GK_points_CumHaz[[i]])
    for(j in 1:upred.lngth) {
      newdata.GK.CumHaz[[i]][[j]] <- right_rows(newdata[newdata[[idVar]] == i, ],
                                                newdata[newdata[[idVar]] == i, ][[timeVar]], 
                                                newdata[newdata[[idVar]] == i, ][[idVar]], 
                                                rbind(GK_points_CumHaz[[i]][j, ]))
      newdata.GK.CumHaz[[i]][[j]][[timeVar]] <- c(t(rbind(GK_points_CumHaz[[i]][j, ])))
      newdata.GK.CumHaz[[i]][[j]][[idVar]] <- match(newdata.GK.CumHaz[[i]][[j]][[idVar]], 
                                                    unique(newdata.GK.CumHaz[[i]][[j]][[idVar]]))
    }
  }
  postMeans <- object$statistics$postMeans
  betas <- postMeans[grep("betas", names(postMeans), fixed = TRUE)] 
  sigmas <- postMeans[grep("sigma", names(postMeans), fixed = TRUE)]
  sigma <- vector("list", n_outcomes)
  if (any(which_gaussian <- which(fams == "gaussian"))) {
    sigma[which_gaussian] <- sigmas
  }
  D <- postMeans[grep("^D", names(postMeans), fixed = FALSE)]
  gammas <- postMeans[grep("^gammas", names(postMeans), fixed = FALSE)]
  alphas <- postMeans[grep("^alphas", names(postMeans), fixed = FALSE)]
  Bs_gammas <- postMeans[grep("^Bs.*gammas$", names(postMeans), fixed = FALSE)]
  invD <- postMeans[grep("inv_D", names(postMeans))]
  Formulas <- object$model_info$Formulas
  #nams_Formulas <- names(Formulas)
  #nams_Formulas <- gsub("_[^_]+$", "", nams_Formulas)
  #outcome <- match(nams_Formulas, respVars)
  #outcome <- sapply(respVars, grep, x = names(Formulas), fixed = TRUE)
  find_outcome <- function (nams_Formulas, repsVars) {
    out <- numeric(length(nams_Formulas))
    for (i in seq_along(respVars)) {
      ind <- grep(respVars[i], nams_Formulas, fixed = TRUE)
      out[ind] <- i
    }
    out
  }
  outcome <- find_outcome(names(Formulas), respVars)
  indFixed <- lapply(Formulas, "[[", "indFixed")
  indRandom <- lapply(Formulas, "[[", "indRandom")
  RE_inds <- object$model_info$RE_inds
  RE_inds2 <- object$model_info$RE_inds2
  Interactions <- object$model_info$Interactions
  trans_Funs <- object$model_info$transFuns
  survMats.last <- vector("list", n.NewL)
  for (i in seq_len(n.NewL)) {
    newdata.GK.postRE.i <- newdata.GK.postRE[newdata.GK.postRE[[idVar]] == i, ]
    newdata.i <- newdata[newdata[[idVar]] == i, ]
    idL.i <- match(newdata.i[[idVar]], unique(newdata.i[[idVar]]))
    idL2 <- rep(list(unique(idL.i)), n_outcomes)
    idL <- rep(list(idL.i), n_outcomes)
    for(j in 1:length(idL)) {
      idL[[j]] <- idL[[j]][na.inds2[[i]][[j]]]
    }
    idGK <- match(newdata.GK.postRE.i[[idVar]], unique(newdata.GK.postRE.i[[idVar]]))
    ids <- rep(list(idGK), n_outcomes)
    idTs <- ids[outcome]
    newdata.i.id <- last_rows(newdata.i, idL.i)
    newdata.i.id[[timeVar]] <- last.time[[i]]
    Us <- lapply(TermsU, function(term) {
      model.matrix(term, data = newdata.GK.postRE.i)
    })
    mfX <- lapply(TermsL[grep("TermsX", names(TermsL))], 
                  FUN = function (x) model.frame.default(delete.response(x), data = newdata.GK.postRE.i))
    mfX_long <- lapply(TermsL[grep("TermsX", names(TermsL))], 
                       FUN = function (x) model.frame.default(delete.response(x), data = newdata.i))
    mfX_long.resp <- lapply(TermsL[grep("TermsX", names(TermsL))], 
                            FUN = function (x) model.frame.default(x, data = newdata.i))
    y_long <- lapply(mfX_long.resp, model.response)
    mfZ <- lapply(TermsL[grep("TermsZ", names(TermsL))], 
                  FUN = function (x) model.frame.default(x, data = newdata.GK.postRE.i))
    mfZ_long <- lapply(TermsL[grep("TermsZ", names(TermsL))], 
                       FUN = function (x) model.frame.default(x, data = newdata.i))
    X_surv_H_postRE <- mapply(FUN = function (x, y) model.matrix.default(x, y), 
                              formulasL2[grep("TermsX", names(formulasL2))], mfX, SIMPLIFY = FALSE)
    X_long <- mapply(FUN = function (x, y) model.matrix.default(x, y), 
                     formulasL2[grep("TermsX", names(formulasL2))], mfX_long.resp, SIMPLIFY = FALSE)
    Xbetas_postRE <- Xbetas_calc(X_long, betas)
    Z_surv_H_postRE <- mapply(FUN = function (x, y) model.matrix.default(x, y), 
                              formulasL2[grep("TermsZ", names(formulasL2))], mfZ, SIMPLIFY = FALSE)
    Z_long <- mapply(FUN = function (x, y) model.matrix.default(x, y), 
                     formulasL2[grep("TermsZ", names(formulasL2))], mfZ_long, SIMPLIFY = FALSE)
    for (j in seq_len(length(Z_long))) {
      Z_long[[j]] <- Z_long[[j]][na.inds2[[i]][[j]], , drop = FALSE]
    }
    XXs <- build_model_matrix(TermsFormulas_fixed, newdata.GK.postRE.i)
    XXsbetas <- Xbetas_calc(XXs, betas, indFixed, outcome)
    ZZs <- build_model_matrix(TermsFormulas_random, newdata.GK.postRE.i)
    W1s <- splines::splineDesign(control$knots, GK_points_postRE[i, ], ord = control$ordSpline, outer.ok = TRUE)
    W2s <- model.matrix(formulasS, newdata.GK.postRE.i)[, -1, drop = FALSE]
    post_b_input_only <- lapply(indRandom, function(x) rbind(rep(0, max(x))))
    Wlongs <- designMatLong(XXs, betas, ZZs, post_b_input_only, ids, outcome,
                            indFixed, indRandom, Us, trans_Funs)
    Pw <- unlist(P[idGK]) * wk
    P <- unlist(P)
    col_inds = attr(Wlongs, "col_inds")
    row_inds_Us = seq_len(nrow(Wlongs))
    survMats.last[[i]][["idL"]] <- idL
    survMats.last[[i]][["idL2"]] <- lapply(idL, length)
    survMats.last[[i]][["idL3"]] <- unlist(idL)
    survMats.last[[i]][["idL3"]] <- c(unlist(idL)[-length(unlist(idL))] != unlist(idL)[-1L], TRUE)
    survMats.last[[i]][["idGK"]] <- idGK
    survMats.last[[i]][["idGK_fast"]] <- c(idGK[-length(idGK)] != idGK[-1L], TRUE)
    survMats.last[[i]][["idTs"]] <- idTs
    survMats.last[[i]][["Us"]] <- Us
    factor2numeric <- function (x) {
      if (is.factor(x)) {
        if (length(levels(x)) > 2){
          stop("Currently only binary outcomes can be considered.")
        }
        as.numeric(x == levels(x)[2L])
      } else x
    }
    survMats.last[[i]][["y_long"]] <- lapply(y_long, factor2numeric)
    survMats.last[[i]][["Xbetas"]] <- Xbetas_postRE
    survMats.last[[i]][["Z_long"]] <- Z_long
    survMats.last[[i]][["X_long"]] <- X_long
    survMats.last[[i]][["RE_inds"]] <- RE_inds
    survMats.last[[i]][["RE_inds2"]] <- RE_inds2
    survMats.last[[i]][["W1s"]] <- W1s
    survMats.last[[i]][["W2s"]] <- W2s
    survMats.last[[i]][["XXs"]] <- XXs
    survMats.last[[i]][["XXsbetas"]] <- XXsbetas
    survMats.last[[i]][["ZZs"]] <- ZZs
    survMats.last[[i]][["col_inds"]] <- col_inds
    survMats.last[[i]][["row_inds_Us"]] <- row_inds_Us
    survMats.last[[i]][["Pw"]] <- Pw
    survMats.last[[i]][["P"]] <- P
  }
  survMats <- vector("list", n.NewL)
  for (i in 1:n.NewL) {
    survMats[[i]] <- vector("list", nrow(GK_points_CumHaz[[i]]))
  }
  for (i in 1:n.NewL) {
    upred.lngth <- nrow(GK_points_CumHaz[[i]])
    for (j in 1:upred.lngth) {
      newdata.GK.CumHaz.ij <- newdata.GK.CumHaz[[i]][[j]]
      newdata.i <- newdata[newdata[[idVar]] == i, ]
      idL.i <- match(newdata.i[[idVar]], unique(newdata.i[[idVar]]))
      idGK <- match(newdata.GK.CumHaz.ij[[idVar]], unique(newdata.GK.CumHaz.ij[[idVar]]))
      ids <- rep(list(idGK), n_outcomes)
      idTs <- ids[outcome]
      newdata.i.id <- last_rows(newdata.i, idL.i)
      newdata.i.id[[timeVar]] <- last.time[[i]]
      Us <- lapply(TermsU, function(term) {
        model.matrix(term, data = newdata.GK.CumHaz.ij)
      })
      XXs <- build_model_matrix(TermsFormulas_fixed, newdata.GK.CumHaz.ij)
      XXsbetas <- Xbetas_calc(XXs, betas, indFixed, outcome)
      ZZs <- build_model_matrix(TermsFormulas_random, newdata.GK.CumHaz.ij)
      W1s <- splines::splineDesign(control$knots, GK_points_CumHaz[[i]][j, ], ord = control$ordSpline, outer.ok = TRUE)
      W2s <- model.matrix(formulasS, newdata.GK.CumHaz.ij)[, -1, drop = FALSE]
      post_b_input_only <- lapply(indRandom, function(x) rbind(rep(0, max(x))))
      Wlongs <- designMatLong(XXs, betas, ZZs, post_b_input_only, ids, outcome,
                              indFixed, indRandom, Us, trans_Funs)
      Pw <- unlist(P2[[i]][[j]][idGK]) * wk
      col_inds = attr(Wlongs, "col_inds")
      row_inds_Us = seq_len(nrow(Wlongs))
      survMats[[i]][[j]][["idGK_fast"]] <- c(idGK[-length(idGK)] != idGK[-1L], TRUE)
      survMats[[i]][[j]][["idTs"]] <- idTs
      survMats[[i]][[j]][["Us"]] <- Us
      survMats[[i]][[j]][["RE_inds2"]] <- RE_inds2
      survMats[[i]][[j]][["W1s"]] <- W1s
      survMats[[i]][[j]][["W2s"]] <- W2s
      survMats[[i]][[j]][["XXs"]] <- XXs
      survMats[[i]][[j]][["XXsbetas"]] <- XXsbetas
      survMats[[i]][[j]][["ZZs"]] <- ZZs
      survMats[[i]][[j]][["col_inds"]] <- col_inds
      survMats[[i]][[j]][["row_inds_Us"]] <- row_inds_Us
      survMats[[i]][[j]][["Pw"]] <- Pw
      survMats[[i]][[j]][["P2"]] <- P2[[i]][j]
    }
  }
  modes.b <- matrix(0, n.NewL, ncol(D[[1]]))
  invVars.b <- Vars.b <- vector("list", n.NewL)
  for (i in 1:n.NewL) {
    Data <- list("idL" = survMats.last[[i]]$idL, "idL2" = survMats.last[[i]]$idL2, 
                 "idGK" = which(survMats.last[[i]]$idGK_fast) - 1, "ids" = survMats.last[[i]]$ids, 
                 "idTs" = survMats.last[[i]]$idTs, "Us" = survMats.last[[i]]$Us, 
                 "y" = survMats.last[[i]]$y_long, "Xbetas" = survMats.last[[i]]$Xbetas, 
                 "Z" = survMats.last[[i]]$Z_long, "RE_inds" = survMats.last[[i]]$RE_inds, 
                 "RE_inds2" = survMats.last[[i]]$RE_inds2, "W1s" = survMats.last[[i]]$W1s, 
                 "W2s" = survMats.last[[i]]$W2s, "XXsbetas" = survMats.last[[i]]$XXsbetas, 
                 "ZZs" = survMats.last[[i]]$ZZs, "col_inds"= survMats.last[[i]]$col_inds, 
                 "row_inds_Us" = survMats.last[[i]]$row_inds_Us, "Bs_gammas" = unlist(Bs_gammas), 
                 "gammas" = unlist(gammas), "alphas" = unlist(alphas), "fams" = fams, 
                 "links" = links, "sigmas" = sigma, "invD" = invD[[1]], 
                 "Pw" = survMats.last[[i]]$Pw, "trans_Funs" = trans_Funs, 
                 "wk" = wk, 
                 "idL3" = which(survMats.last[[i]]$idGK_fast) - 1)
    if (is.null(Data$gammas))
      Data$gammas <- numeric(0)
    ff <- function (b, Data) -log_post_RE_svft(b, Data = Data)
    gg <- function (b, Data) cd(b, ff, Data = Data, eps = 1e-03)
    
    t_start = Sys.time()
    start <- rep(0, ncol(D[[1]]))
    opt <- optim(start, ff, gg, Data = Data, method = "BFGS", hessian = TRUE, 
                 control = list(maxit = 200, parscale = rep(0.1, ncol(D[[1]]))))
    t_end = Sys.time()
    
    modes.b[i, ] <- opt$par
    
    invVars.b[[i]] <- opt$hessian / scale
    Vars.b[[i]] <- scale * solve(opt$hessian)
  }
  set.seed(seed)
  out <- vector("list", M)
  success.rate <- matrix(FALSE, M, n.NewL)
  b.old <- b.new <- modes.b
  mcmc <- object$mcmc
  mcmc <- mcmc[names(mcmc) != "b"]
  if (M > nrow(mcmc$betas1)) {
    warning("'M' cannot be set greater than ", nrow(mcmc$betas1))
    M <- nrow(mcmc$betas1)
    out <- vector("list", M)
    success.rate <- matrix(FALSE, M, n.NewL)
  }
  samples <- sample(nrow(mcmc$betas1), M)
  mcmc[] <- lapply(mcmc, function (x) {
    if (length(dim(x)) == 3) {
      x[samples, , ]
    } else if (length(dim(x)) == 2) {
      x[samples, , drop = FALSE]
    } else if (all(is.null(dim(x)), length(x) > 0)) {
      x[samples]
    } else {
      x[]
    }
  })
  proposed.b <- mapply(rmvt, mu = split(modes.b, row(modes.b)), Sigma = Vars.b, 
                       MoreArgs = list(n = M, df = 4), SIMPLIFY = FALSE)
  proposed.b[] <- lapply(proposed.b, function (x) if (is.matrix(x)) x else rbind(x))
  dmvt.proposed <- mapply(dmvt, x = proposed.b, mu = split(modes.b, row(modes.b)), 
                          Sigma = Vars.b, MoreArgs = list(df = 4, log = TRUE), 
                          SIMPLIFY = FALSE)
  SS <- vector("list", n.NewL)
  for (i in 1:n.NewL) {
    for (m in 1:M) {
      betas.new <- lapply(mcmc[grep("betas", names(mcmc))], function (x) x[m, ])
      XXsbetas.new <- Xbetas_calc(survMats.last[[i]][["XXs"]], betas.new, indFixed, outcome)
      alphas.new <- mcmc$alphas[m, ]
      invD.new <- mcmc$inv_D[m, ,] 
      Bs_gammas.new <- mcmc$Bs_gammas[m, ]
      gammas.new <- mcmc$gammas[m, ]
      sigma <- lapply(mcmc[grep("sigma", names(mcmc))], function (x) x[m])
      sigma.new <- vector("list", n_outcomes)
      if (any(which_gaussian <- which(fams == "gaussian"))) {
        sigma.new[which_gaussian] <- sigma
      }
      p.b <- proposed.b[[i]][m, ]
      dmvt.old <- dmvt(b.old[i, ], modes.b[i, ], invSigma = invVars.b[[i]], df = 4, log = TRUE)
      dmvt.prop <- dmvt.proposed[[i]][m]
      Xbetasnew.m <- Xbetas_calc(survMats.last[[i]]$X_long, betas.new)
      Data <- list("idL" = survMats.last[[i]]$idL, "idL2" = survMats.last[[i]]$idL2, 
                   "idGK" = which(survMats.last[[i]]$idGK_fast) - 1, "ids" = survMats.last[[i]]$ids, 
                   "idTs" = survMats.last[[i]]$idTs, "Us" = survMats.last[[i]]$Us, 
                   "y" = survMats.last[[i]]$y_long, "Xbetas" = Xbetasnew.m, 
                   "Z" = survMats.last[[i]]$Z_long, "RE_inds" = survMats.last[[i]]$RE_inds, 
                   "RE_inds2" = survMats.last[[i]]$RE_inds2, "W1s" = survMats.last[[i]]$W1s, 
                   "W2s" = survMats.last[[i]]$W2s, "XXsbetas" = XXsbetas.new, 
                   "ZZs" = survMats.last[[i]]$ZZs, "col_inds"= survMats.last[[i]]$col_inds, 
                   "row_inds_Us" = survMats.last[[i]]$row_inds_Us, "Bs_gammas" = Bs_gammas.new, 
                   "gammas" = gammas.new, "alphas" = alphas.new, "fams" = fams, 
                   "links" = links, "sigmas" = sigma.new, "invD" = invD.new, 
                   "Pw" = survMats.last[[i]]$Pw, "trans_Funs" = trans_Funs, 
                   "wk" = wk, "idL3" = which(survMats.last[[i]]$idGK_fast) - 1)
      if (is.null(Data$gammas))
        Data$gammas <- numeric(0)
      a <- min(exp(log_post_RE_svft(p.b, Data = Data) + dmvt.old - 
                     log_post_RE_svft(b.old[i, ], Data = Data) - dmvt.prop), 1)
      ind <- runif(1) <= a
      if (!is.na(ind) && ind) {
        b.new[i, ] <- p.b
      }
      logS.pred <- numeric(length(times.to.pred_upper[[i]]))
      for (l in seq_along(logS.pred)) {
        XXsbetas.newl <- Xbetas_calc(survMats[[i]][[l]][["XXs"]], betas.new, indFixed, outcome)
        Datal <- list("idGK" = which(survMats[[i]][[l]][["idGK_fast"]]) - 1, "idTs" = survMats[[i]][[l]][["idTs"]], 
                      "Us" = survMats[[i]][[l]][["Us"]], "RE_inds2" = survMats[[i]][[l]][["RE_inds2"]], 
                      "W1s" = survMats[[i]][[l]][["W1s"]], "W2s" = survMats[[i]][[l]][["W2s"]], 
                      "XXsbetas" = XXsbetas.newl, "ZZs" = survMats[[i]][[l]][["ZZs"]], 
                      "col_inds" = survMats[[i]][[l]][["col_inds"]], 
                      "row_inds_Us" = survMats[[i]][[l]][["row_inds_Us"]], 
                      "Bs_gammas" = Bs_gammas.new, 
                      "gammas" = gammas.new, "alphas" = alphas.new, 
                      "Pw" = survMats[[i]][[l]][["Pw"]], "trans_Funs" = trans_Funs, 
                      "wk" = wk)
        if (is.null(Datal$gammas))
          Datal$gammas <- numeric(0)
        logS.pred[l] <- survPred_svft_2(b.new[i, ], Data = Datal)
      }
      SS[[i]] <- if (log) logS.pred else exp(cumsum(logS.pred))
      b.old <- b.new
      out[[m]][[i]] <- SS[[i]]
    }
  }
  res <- vector("list", n.NewL)
  for (i in seq_len(n.NewL)) {
    rr <- sapply(out, "[[", i)
    if (!is.matrix(rr)) {
      rr <- rbind(rr)
    }
    res[[i]] <- cbind(
      times = times.to.pred_upper[[i]], 
      "Mean" = rowMeans(rr, na.rm = TRUE), 
      "Median" = apply(rr, 1L, median, na.rm = TRUE), 
      "Lower" = apply(rr, 1L, quantile, probs = CI.levels[1], na.rm = TRUE),
      "Upper" = apply(rr, 1L, quantile, probs = CI.levels[2], na.rm = TRUE)
    )
    rownames(res[[i]]) <- as.character(seq_len(NROW(res[[i]])))
  }
  y. <- lapply(y, split, id)
  y <- vector("list", n.NewL)
  for (i in seq_len(n.NewL)) {
    y[[i]] <- lapply(y., "[[", i)
  }
  newdata. <- do.call(rbind, mapply(function (d, t) {
    if (d[[timeVar]][nrow(d)] < t) {
      d. <- rbind(d, d[nrow(d), ])
      d.[[timeVar]][nrow(d.)] <- t
      d.
    } else d
  }, split(newdata, id.), last.time, SIMPLIFY = FALSE))
  id. <- newdata.[[idVar]]
  id. <- match(id., unique(id.))
  mfX. <- lapply(TermsL[grep("TermsX", names(TermsL))], 
                 FUN = function (x) model.frame.default(delete.response(x), data = newdata.))
  mfZ. <- lapply(TermsL[grep("TermsZ", names(TermsL))], 
                 FUN = function (x) model.frame.default(x, data = newdata.))
  X. <- mapply(FUN = function (x, y) model.matrix.default(x, y), 
               formulasL2[grep("TermsX", names(formulasL2))], mfX., SIMPLIFY = FALSE)
  Z. <- mapply(FUN = function (x, y) model.matrix.default(x, y), 
               formulasL2[grep("TermsZ", names(formulasL2))], mfZ., SIMPLIFY = FALSE)
  fitted.y <- vector("list", length(X.))
  for (i in seq_len(n.NewL)) {
    fits <- vector("list", length(X.))
    for (j in seq_along(X.)) {
      id_i <- id. == i
      fits[[j]] <- c(X.[[j]][id_i, , drop = FALSE] %*% betas[[j]]) + 
        rowSums(Z.[[j]][id_i, , drop = FALSE] * modes.b[rep(i, sum(id_i)), RE_inds[[j]], drop = FALSE])
    }
    fitted.y[[i]] <- fits
  }
  fitted.y <- fitted.y[!sapply(fitted.y, is.null)]
  fitted.times <- split(newdata.[[timeVar]], factor(newdata.[[idVar]]))
  names(res) <- names(last.time) <- names(obs.times) <- names(fitted.times) <- names(y) <- names(fitted.y) <- levels(idNewL)
  y[] <- lapply(y, function (x, nam) {names(x) <- nam; x}, nam = respVars)
  fitted.y[] <- lapply(fitted.y, function (x, nam) {names(x) <- nam; x}, nam = respVars)
  names(families) <- respVars
  res <- list(summaries = res, full.results = out, survTimes = survTimes, last.time = last.time, 
              obs.times = obs.times, y = y, M = M, families = families, respVars = respVars,
              fitted.times = fitted.times, 
              fitted.y = fitted.y, ry = lapply(componentsL[grep("y", names(componentsL))], FUN = range, na.rm = TRUE), 
              nameYs = unname(lapply(formulasL[grep("TermsX", names(formulasL))], FUN = function(x) paste(x)[2L])))
  rm(list = ".Random.seed", envir = globalenv())
  class(res) <- "survfit.mvJMbayes"
  res
  
  return(opt$par)
}

environment(survfit_temp) <- asNamespace('JMbayes')
