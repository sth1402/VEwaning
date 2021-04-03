#' @include ulbound.R
pre_dNYb <- function(data, times, KRfR, fEX, lag, dObj, delta, minWgt, maxWgt) {

  # I(E_i < t_m) I(U_i >= t_m) i = 1:n, m = 1:M
  Yt <- outer(X = data[,dObj$E], Y = times, FUN = "<") *
        outer(X = data[,dObj$U], Y = {times-1e-8}, FUN = ">=")

  # I(U_i == t_m) Delta_i
  dNt <- outer(X = data[,dObj$U], Y = {times-1e-8}, FUN = ">=") *
         outer(X = data[,dObj$U], Y = {times+1e-8}, FUN = "<=") *
         delta

  #  Construct numerator indicators for dtilNb and tilYb

  irt <- outer(X = data[,dObj$R], Y = {times-1e-8}, FUN = ">=")

  # (1-A_i)I(R_i >= t_m)
  A1IRt <- {1L - data[,dObj$A]} * irt
           

  # A_i I(E_i + lag <= t_m) I(R_i >= t_m)
  AIEelltR <- data[,dObj$A] * irt *
              outer(X = data[,dObj$E] + lag, Y = {times+1e-8}, FUN = "<=")

  #  stabilized weights

  if (is.null(x = fEX)) {
    w0h0 <- 1.0
    w1h1 <- 1.0
  } else {
    w0h0 <- ulbound(wgt = fEX * KRfR$KR.0.stab, 
                    minWgt = minWgt,  
                    maxWgt = maxWgt)
    w1h1 <- ulbound(wgt = fEX * KRfR$KR.1.stab, 
                    minWgt = minWgt,  
                    maxWgt = maxWgt)
  }
    
  ##  dtilNb(t) and Ytilb(t)

  # (1-A_i)I(R_i >= t_m) fEX KRfR0
  wt0 <- A1IRt * w0h0

  # A_i I(E_i + lag <= t_m) I(R_i >= t_m) fEx KRfR1
  wt1 <- AIEelltR * w1h1
    
  dNtTilde <- dNt * {wt0 + wt1}

  return( list("dNtTilde" = dNtTilde, "wt0" = wt0, "wt1" = wt1, "Y" = Yt) )
}
