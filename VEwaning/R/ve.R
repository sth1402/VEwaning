#' Retrieve the Estimated Vaccine Efficacy
#'
#' Uses a prior veWaning() analysis to estimate the vaccine efficacy
#'   at the provided times since vaccination.
#'
#' @param x An object of class VEwaning. The object returned by a call to
#'   veWaning()
#'
#' @param taus A numeric vector object. The times since vaccination at which
#'   the vaccine efficacy is to be estimated.
#'
#' @returns A matrix object. The first column contains the times since
#'   vaccination at which the estimates are provided; the second column
#'   contains estimated vaccine efficacy; and the third is the standard error.
#'
#'
#' @name ve
#' @examples
#' data(veWaningData)
#'
#' set.seed(1234)
#'
#' ind <- sample(1:nrow(veWaningData), 2500)
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
#'                 v = c(20))
#'
#' ve(x = res, taus = c(10,20,30,40,50))
#' @export 
ve <- function(x, taus) {

  lag <- attr(x = x, which = "lag")
  maxTau <- attr(x = x, which = "maxTau")
  gFunc <- attr(x = x, which = "gFunc")
  v <- attr(x = x, which = "v")
  if (is.null(x = v)) v <- 1.0

  if (any(taus < lag)) {
    message("tau values < lag have been removed")
    taus <- taus[taus >= lag]
    if (length(x = taus) == 0L) {
      stop("inappropriate tau values provided", call. = FALSE)
    }
  }

  if (any(taus > maxTau)) {
    message("tau values > than the maximum value in analysis data are ignored")
    taus <- taus[taus <= maxTau]
    if (length(x = taus) == 0L) {
      stop("inappropriate tau values provided", call. = FALSE)
    }
  }

  times <- matrix(data = taus - lag, ncol = 1L)
  gFuncR <- gFunction(gFunc = gFunc, u = times, theta = x$theta, knts = v)
  rate <- exp(x = x$theta[1L] + gFuncR$gu[,1L])

  drate <- matrix(data = 0.0, 
                  nrow = length(x = taus), 
                  ncol = length(x = gFuncR$gutheta))

  for (i in 1L:ncol(x = drate)) {
    drate[,i] <- gFuncR$gutheta[[ i ]]*rate
  }

  se <- diag(x = drate %*% x$cov %*% t(x = drate))

  return( cbind("tau" = taus, "VE" = 1.0 - rate, "SE" = sqrt(x = se)) )

}
