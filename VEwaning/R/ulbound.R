ulbound <- function(wgt, minWgt, maxWgt) {

  wgt[wgt > maxWgt] <- maxWgt
  wgt[wgt < minWgt] <- minWgt

  return( wgt )

}
