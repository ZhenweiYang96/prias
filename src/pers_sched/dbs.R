dbs_deg1 = function (x, df = NULL, knots = NULL, intercept = FALSE, Boundary.knots = range(x), 
          eps = 0.001) 
{
  ns.x <- if (is.null(knots)) {
    splines::bs(x, df = df, degree=1, intercept = intercept, Boundary.knots = Boundary.knots)
  }
  else {
    splines::bs(x, knots = knots, degree=1, intercept = intercept, 
                Boundary.knots = Boundary.knots)
  }
  kn <- attr(ns.x, "knots")
  Bkn <- attr(ns.x, "Boundary.knots")
  ex <- pmax(abs(x), 1)
  x1 <- x + eps * ex
  x2 <- x - eps * ex
  ns.xeps1 <- splines::bs(x1, knots = kn, Boundary.knots = Bkn, degree=1,
                          intercept = intercept)
  ns.xeps2 <- splines::bs(x2, knots = kn, Boundary.knots = Bkn, degree=1,
                          intercept = intercept)
  out <- (ns.xeps1 - ns.xeps2)/c(x1 - x2)
  attr(out, "eps") <- eps
  attr(out, "class") <- c("dbs", "basis", "matrix")
  out
}