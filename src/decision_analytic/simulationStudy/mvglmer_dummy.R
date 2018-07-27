mvglmer_dummy <- function (formulas, data, families, engine = c("JAGS", "STAN"), 
                     overdispersion = FALSE, priors = NULL, init = NULL, 
                     control = NULL, ...) {
  cl <- match.call()
  engine <- match.arg(engine)
  if (!is.list(families))
    stop("'families' must be a list of family objects.")
  # depending on the input of the user, set families to the corresponding
  # base R functions
  families[] <- lapply(families, function (family) {
    if (is.character(family))
      family <- get(family, mode = "function", envir = parent.frame())
    if (is.function(family))
      family <- family()
    family
  })
  # using the formulas and the data extract and constuct the objects required to
  # pass to JAGS
  components <- lapply(unname(formulas), JMbayes:::extractFrames, data = data)
  # perform some checks
  if (!all((ns <- sapply(components, `[[`, 'n')) == components[[1L]][['n']])) {
    stop("it seems that the number of subjects differ between the responses. More ",
         "specifically, according to the data the number of subjects per response is ",
         paste(paste(sapply(components, `[[`, 'respVar'), ns, sep = " = "),
               collapse = ", "), ".\n")
  }
  components <- unlist(components, recursive = FALSE)
  n_outcomes <- length(formulas)
  names(components) <- paste0(names(components),
                              rep(seq_len(n_outcomes),
                                  each = length(components) / n_outcomes))
  colmns_HC <- components[grep("colmns_HC", names(components), fixed = TRUE)]
  colmns_nHC <- components[grep("colmns_nHC", names(components), fixed = TRUE)]
  seq_outcomes <- seq_len(n_outcomes)
  nams_vars <- c("N", "id", "Z", "Xhc", "ncx", "y")
  vars <- paste0(rep(nams_vars, each = n_outcomes), seq_outcomes)
  if (any(ind_td <- sapply(colmns_nHC, length))) {
    vars <- c(vars, paste0("X", which(ind_td > 0)))
  }
  Data <- c(list(n = components$n1), components[vars])
  Data$n_RE <- sum(unlist(components[grep("ncz", names(components), fixed = TRUE)]))
  RE_inds <- mapply(function (sq, incr) seq_len(sq) + incr,
                    sq = components[grep("ncz", names(components), fixed = TRUE)],
                    incr = cumsum(c(0, head(sapply(colmns_HC, length), -1))),
                    SIMPLIFY = FALSE)
  names(RE_inds) <- paste0("RE_ind", seq_along(RE_inds))
  Data <- c(Data, RE_inds, unlist(colmns_HC, recursive = FALSE), colmns_nHC)
  # control
  con <- list(n.processors = parallel::detectCores() - 1, n.chains = 2,
              working.directory = getwd(), clear.model = TRUE,
              seed = 1L, optimize_only = FALSE, verbose = FALSE)
  if (engine == "JAGS") {
    con$n.iter <- 28000L
    con$n.burnin <- 3000L
    con$n.thin <- 50L
    con$n.adapt <- 3000L 
  } else {
    con$n.iter <- 1000
    con$n.warmup <- floor(con$n.iter / 2)
    con$n.thin <- 1
    con$adapt_delta <- 0.8
  }
  control <- c(control, list(...))
  namC <- names(con)
  con[(namc <- names(control))] <- control
  if (!any(namc == "n.thin")) {
    con$n.thin <- if (engine == "JAGS") {
      max(1, floor((con$n.iter - con$n.burnin) * con$n.chains / 1000))
    } else {
      max(1, floor((con$n.iter - con$n.warmup) * con$n.chains / 1000))
    }
  }
  if (length(noNms <- namc[!namc %in% namC]) > 0)
    warning("unknown names in control: ", paste(noNms, collapse = ", "))
  ######################################################################################
  # Priors
  if (engine == "JAGS") {
    prs <- list(priorR_D = diag(rep(as.numeric(NA), Data$n_RE), Data$n_RE),
                priorK_D = Data$n_RE + 1, A_RD = 0.5, B_RD = 0.01,
                tau_half_cauchy = 0.1)
    pr_taus_betas <- rep(list(0.01), n_outcomes)
    names(pr_taus_betas) <- paste0("tau_betas", seq_len(n_outcomes))
    prs <- c(prs, pr_taus_betas)
    if (any(sapply(families, function (x) x$family == "gaussian")) || overdispersion) {
      prs$A_tau <- 0.01
      prs$B_tau <- 0.01
    }
  } else {
    prs <- list(scale_sigmas = 5, scale_diag_D = 3, lkj_shape = 2,
                priorK_D = Data$n_RE + 1)
    pr_scale_betas <- rep(list(10), n_outcomes)
    names(pr_scale_betas) <- paste0("scale_betas", seq_len(n_outcomes))
    prs <- c(prs, pr_scale_betas)
  }
  if (!is.null(priors)) {
    lngths <- lapply(prs[(nam.prs <- names(priors))], length)
    if (!is.list(priors) || !isTRUE(all.equal(lngths, lapply(priors, length)))) {
      warning("'priors' is not a list with elements numeric vectors of appropriate ",
              "length; default priors are used instead.\n")
    } else {
      prs[nam.prs] <- priors
    }
  }
  Data <- c(Data, prs)
  ######################################################################################
  
  #closeAllConnections()
  # parameters to save
  params <- paste0('betas', seq_len(n_outcomes))
  if (any(ind_gs <- sapply(families, function (x) x$family == "gaussian"))) {
    params <- c(params, paste0("sigma", which(ind_gs)))
  }
  params <- if (engine == "JAGS") c(params, "inv_D", "b") else c(params, "D", "b")
  inits <- function () {
    ints <- lapply(components[grep("ncx", names(components), fixed = TRUE)],
                   rnorm, sd = 0.1)
    names(ints) <- paste0('betas', seq_len(n_outcomes))
    ints$u <- drop(matrix(rnorm(Data$n * Data$n_RE), Data$n,
                          Data$n_RE))
    ints$inv_D <- ints$D <- if (Data$n_RE > 1) diag(Data$n_RE) else 1
    if (any(ind_gs)) {
      nms <- which(ind_gs)
      taus <- rep(list(1), length(nms))
      names(taus) <- paste0("tau", nms)
      ints <- c(ints, taus)
    }
    ints
  }
  
  out <- list(components = components, data = data)

  class(out) <- "mvglmer"
  out
}

