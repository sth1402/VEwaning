#' Print Analysis Results
#'
#' Print the primary results of the analysis
#'
#' @param x An object of class VEwaning. The object returned by a call to
#'   veWaning()
#'
#' @param ... Ignored
#'
#'
#' @name print
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
#'                 v = c(15,30))
#'
#' print(x = res)
#' @method print VEwaning
#' @export 
print.VEwaning <- function(x, ...) {

  attr(x, "gFunc") <- NULL
  attr(x, "maxTau") <- NULL
  attr(x, "lag") <- NULL
  attr(x, 'v') <- NULL

  x <- unclass(x = x)

  print(x)

}
