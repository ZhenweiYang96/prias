#Adding 2 new functions
fixef_weighted <- function (object){
  if (!inherits(object, "mvJMbayes"))
    stop("Use only with 'mvJMbayes' objects.\n")
  comps <- object$model_info$mvglmer_components
  nams_outcomes <- unlist(comps[grep("respVar", names(comps), fixed = TRUE)], 
                          use.names = FALSE)
  pMeans <- object$statistics$postwMeans
  betas <- pMeans[grep("betas", names(pMeans), fixed = TRUE)]
  names(betas) <- nams_outcomes
  betas
}

ranef_weighted <- function (object, as_list=FALSE){
  if (!inherits(object, "mvJMbayes"))
    stop("Use only with 'mvJMbayes' objects.\n")
  out <- object$statistics$postwMeans$b
  dimnames(out) <- list(names(object$model_info$coxph_components$Time), 
                        colnames(object$statistics$postwMeans$D))
  if (as_list) {
    components <- object$model_info$mvglmer_components
    ncols <- components[grep("ncz", names(components), fixed = TRUE)]
    RE_inds <- object$model_info$RE_inds
    out <- lapply(RE_inds, function (i) out[, i, drop = FALSE])
  }
  out
}

fitted_weighted <- function (object) {
  if (!inherits(object, "mvJMbayes")) 
    stop("Use only with 'mvJMbayes' objects.\n")
  components <- object$model_info$mvglmer_components
  families <- object$model_info$families
  n_outcomes <- length(families)
  seq_n_outcomes <- seq_len(n_outcomes)
  X <- components[paste0("X", seq_n_outcomes)]
  fitY <- mapply("%*%", X, fixef_weighted(object), SIMPLIFY = FALSE)
  names(fitY) <- row.names(object$Data$data)
  ids <- components[paste0("id", seq_n_outcomes)]
  Zs <- components[paste0("Z", seq_n_outcomes)]
  bs <- ranef_weighted(object, as_list = TRUE)
  Zb_fun <- function (Z, b, id) rowSums(Z * b[id, , drop = FALSE])
  Zbs <- mapply(Zb_fun, Zs, bs, ids, SIMPLIFY = FALSE)
  fitY <- mapply("+", fitY, Zbs, SIMPLIFY = FALSE)
  fitY
}