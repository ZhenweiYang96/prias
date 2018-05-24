qqtdist = function (y, ylim, mu=0, sigma=1, df, main = "T Dist Q-Q Plot", xlab = "Theoretical Quantiles", 
          ylab = "Sample Quantiles", plot.it = TRUE, datax = FALSE, 
          ...) 
{
  if (has.na <- any(ina <- is.na(y))) {
    yN <- y
    y <- y[!ina]
  }
  if (0 == (n <- length(y))) 
    stop("y is empty or has only NAs")
  if (plot.it && missing(ylim)) 
    ylim <- range(y)
  
  
  x <- qt(ppoints(n), df=df)[order(order(y))]
  #x <- JMbayes:::qgt(ppoints(n), mu=mu, sigma=sigma, df=df)[order(order(y))]
  #x <- qnorm(ppoints(n))[order(order(y))]
  if (has.na) {
    y <- x
    x <- yN
    x[!ina] <- y
    y <- yN
  }
  if (plot.it) 
    if (datax) 
      plot(y, x, main = main, xlab = ylab, ylab = xlab, 
           xlim = ylim, ...)
  else plot(x, y, main = main, xlab = xlab, ylab = ylab, 
            ylim = ylim, ...)
  invisible(if (datax) list(x = y, y = x) else list(x = x, 
                                                    y = y))
}

qqlineTdist=function (y, datax = FALSE, distribution = qt, 
                      mu=0, sigma=1, df, probs = c(0.25, 0.75), qtype = 7, ...) 
{
  stopifnot(length(probs) == 2, is.function(distribution))
  y <- quantile(y, probs, names = FALSE, type = qtype, na.rm = TRUE)
  #x <- distribution(probs, mu=mu, sigma=sigma, df=df)
  x <- distribution(probs, df=df)
  if (datax) {
    slope <- diff(x)/diff(y)
    int <- x[1L] - slope * y[1L]
  }
  else {
    slope <- diff(y)/diff(x)
    int <- y[1L] - slope * x[1L]
  }
  
  print(paste0("Intercept = ", int, ", Slope = ", slope))
  
  abline(int, slope, ...)
}