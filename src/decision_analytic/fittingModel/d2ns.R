d2ns <- function (x, df = NULL, knots = NULL, intercept = FALSE, Boundary.knots = range(x), eps = 1e-03) {
  ns.x <- if (is.null(knots)) {
    splines::ns(x, df = df, intercept = intercept, Boundary.knots = Boundary.knots)
  } else {
    splines::ns(x, knots = knots, intercept = intercept, Boundary.knots = Boundary.knots)
  } 
  kn <- attr(ns.x, "knots")
  Bkn <- attr(ns.x, "Boundary.knots")
  ex <- pmax(abs(x), 1)
  h = eps * ex
  x1 <- x + h
  x2 <- x - h
  ns.xeps1 <- splines::ns(x1, knots = kn, Boundary.knots = Bkn, intercept = intercept)
  ns.xeps2 <- splines::ns(x2, knots = kn, Boundary.knots = Bkn, intercept = intercept)
  out <- (ns.xeps1 - 2*ns.x + ns.xeps2) / h^2
  attr(out, "eps") <- eps
  attr(out, "class") <- c("dns", "basis", "matrix")
  out
}

environment(d2ns) <- asNamespace('JMbayes')