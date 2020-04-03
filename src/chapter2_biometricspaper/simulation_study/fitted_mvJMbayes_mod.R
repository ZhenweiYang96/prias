fitted.mvJMbayes_mod <- function (object, process = c("Longitudinal", "longitudinal", "Event", 
                                                  "event"), 
                              type = c("Marginal", "marginal", "Subject", "subject"), ...) {
  if (!inherits(object, "mvJMbayes")) 
    stop("Use only with 'mvJMbayes' objects.\n")
  process <- match.arg(process)
  type <- match.arg(type)
  if (process == "Longitudinal" || process == "longitudinal") {
    components <- object$model_info$mvglmer_components
    families <- object$model_info$families
    n_outcomes <- length(families)
    seq_n_outcomes <- seq_len(n_outcomes)
    X <- components[paste0("X", seq_n_outcomes)]
    fitY <- mapply("%*%", X, fixef(object), SIMPLIFY = FALSE)
    names(fitY) <- row.names(object$Data$data)
    if (type == "Subject" || type == "subject") {
      ids <- components[paste0("id", seq_n_outcomes)]
      Zs <- components[paste0("Z", seq_n_outcomes)]
      random_b <- list(apply(object$mcmc$b, c(1,2), mean))
      
      Zb_fun <- function (Z, b, id){
        rowSums(Z * b[id, , drop = FALSE])
      }
      
      Zbs <- mapply(Zb_fun, Zs, random_b, ids, SIMPLIFY = FALSE)
      fitY <- mapply("+", fitY, Zbs, SIMPLIFY = FALSE)
    }
    fitY
  }
}

environment(fitted.mvJMbayes_mod) <- asNamespace('JMbayes')
