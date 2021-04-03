#' Plot Analysis Results
#'
#' Plot the Estimated Vaccine Efficacy
#'
#' @param x An object of class VEwaning. The object returned by a call to
#'   veWaning()
#'
#' @param y Ignored
#'
#' @param ... Ignored
#'
#' @param xlim A numeric vector object. The minimum and maximum tau values
#'   to include in the plot.
#'
#'
#' @name plot
#' @examples
#' data(veWaningData)
#'
#' set.seed(1234)
#'
#' ind <- sample(1:nrow(x = veWaningData), 2500, FALSE)
#'
#' # NOTE: This sample size is chosen for example only -- larger data sets
#' # should be used.
#'
#' res <- veWaning(data = veWaningData[ind,], 
#'                 L = 52,  
#'                 lag = 6,  
#'                 modelGam1 = ~ X1+X2+A+A:X1+A:X2, 
#'                 modelGam2 = ~ X1+X2, 
#'                 modelEntry = ~ X1+X2, 
#'                 modelPsiGam1 = ~ X1+X2, 
#'                 modelPsiGam2 = ~ X1+X2, 
#'                 gFunc = 'piece', 
#'                 v = c(15,30))
#'
#' plot(x = res)
#' @method plot VEwaning
#' @export 
#' @importFrom graphics plot
#' @import ggplot2
plot.VEwaning <- function(x, y, ..., xlim) {

  lag <- attr(x = x, which = "lag")
  maxTau <- attr(x = x, which = "maxTau")
  gFunc <- attr(x = x, which = "gFunc")
  v <- attr(x = x, which = "v")
  if (is.null(x = v)) v <- 1.0

  if (!missing(x = xlim)) {
    xlim <- sort(x = xlim)
    minTau <- xlim[1]
    maxTau <- xlim[2]
  } else {
    minTau = attr(x = x, which = "lag")
  }

  taus = seq(from = minTau, to = maxTau, length.out = 20)

  times <- matrix(data = taus - lag, ncol = 1L)
  gFuncR <- gFunction(gFunc = gFunc, u = times, theta = x$theta, knts = v)

  dg <- NULL
  for (i in 1L:length(x = gFuncR$gutheta)) {
    dg <- cbind(dg, gFuncR$gutheta[[ i ]])
  }

  se <- sqrt(x = diag(x = dg %*% x$cov %*% t(x = dg)))

  veLow <- 1.0 - exp(x = x$theta[1L] + gFuncR$gu[,1L] - 1.96*se)
  ve <- 1.0 - exp(x = x$theta[1L] + gFuncR$gu[,1L])
  veHigh <- 1.0 - exp(x = x$theta[1L] + gFuncR$gu[,1L] + 1.96*se)

  df <- data.frame(taus, veLow, ve, veHigh)

  ggplot(df, aes(x = taus, y = ve)) +
    geom_path(lwd = 1) +
    geom_ribbon(aes(ymin = veLow, ymax = veHigh), alpha = 0.2, fill = "grey") +
    labs(y = expression(paste("VE(", tau, ")")), x = expression(tau)) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"))

}
