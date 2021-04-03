#' @include ulbound.R
pre_dNYu <- function(data, 
                     times,  
                     KRfR,  
                     fEX,  
                     Psiprob,  
                     lag,  
                     dObj,  
                     delta,
                     minWgt,
                     maxWgt) {

  nTimes <- length(x = times)

  # Observed at-risk indicator at time t_m
  # Y = I(E_i < t_m <= U_i)
  # {n x M}
  Y <- outer(X = data[,dObj$E], Y = times, FUN = "<") *
       outer(X = data[,dObj$U], Y = {times-1e-8}, FUN = ">=")

  # Indicator that a participant is observed to be infected at time t_m
  # dNt = I(U_i == t_m, Delta == 1)
  # {n x M}
  dNt <- outer(X = data[,dObj$U], Y = {times-1e-8}, FUN = ">=") *
         outer(X = data[,dObj$U], Y = {times+1e-8}, FUN = "<=") *
         delta

  # Construct numerator indicators for dtilNb and tilYb
  #
  # Ytilde = Y[(1-A)I(t_m - R_i >= lag)
  #   {I(Gamma = 1, Psi = 1) + I(Gamma = 2, Psi = 1)} w0/h01 exp(g(R+l)) +
  #   A I(t_m > R_i){I(Gamma = 1) + I(Gamma = 2)} w1/h11 exp(g(E+l))]

  # (1-A_i) I(t_m - R >= lag)
  A1ItRl <- {1L - data[,dObj$A]} * 
            outer(X = -data[,dObj$R], Y = lag - times, FUN = ">=")

  # A_i I(t_m > R)
  AItR <- data[,dObj$A]*
          outer(X = data[,dObj$R], Y = times, FUN = "<=")

  # I(Gamma = 1)
  Gam1 <- data[,dObj$Gam] == 1L
  # I(Gamma = 2)
  Gam2 <- data[,dObj$Gam] == 2L

  # I(Gamma = 1, Psi = 1)
  Gam1Psi1 <- Gam1 * {data[,dObj$Psi] == 1L}

  # I(Gamma = 2, Psi = 1)
  Gam2Psi1 <- Gam2 * {data[,dObj$Psi] == 1L}
  
  if (is.null(x = fEX)) {
    # w0/h01 for Gamma = 1
    w0h01Gam1 <- 1.0

    # w0/h01 for Gamma  = 2
    w0h01Gam2 <- 1.0

    # w1/h11 for Gamma = 1
    w1h11Gam1 <- 1.0

    # w1/h11 for Gamma = 2
    w1h11Gam2 <- 1.0
  } else {
    # w0/h01 for Gamma = 1
    w0h01Gam1 <- ulbound(wgt = KRfR$fR1.0.stab * Psiprob$pPsi1.stab * fEX, 
                         minWgt = minWgt,  
                         maxWgt = maxWgt)

    # w0/h01 for Gamma  = 2
    w0h01Gam2 <- ulbound(wgt = KRfR$fR2.0.stab * Psiprob$pPsi2.stab * fEX, 
                         minWgt = minWgt,  
                         maxWgt = maxWgt)

    # w1/h11 for Gamma = 1
    w1h11Gam1 <- ulbound(wgt = KRfR$fR1.1.stab*fEX, 
                         minWgt = minWgt,  
                         maxWgt = maxWgt)

    # w1/h11 for Gamma = 2
    w1h11Gam2 <- ulbound(wgt = KRfR$fR2.1.stab*fEX, 
                         minWgt = minWgt,  
                         maxWgt = maxWgt)
  }

  # I(Gamma = 1, Psi = 1) w0/h01 + I(Gamma = 2, Psi = 1) * w0/h01
  w0h12 <- Gam1Psi1 * w0h01Gam1 + Gam2Psi1 * w0h01Gam2

  # {I(Gamma = 1) w1/h11 + I(Gamma = 2) w1/h11
  w1h12 <- Gam1 * w1h11Gam1 + Gam2 * w1h11Gam2

  # (1-A_i) I(t_m - R >= lag) 
  #   {I(Gamma = 1, Psi = 1) w0/h01 + I(Gamma = 2, Psi = 1) w0/h01}
  wt0 <- w0h12*A1ItRl

  # A_i I(t_m > R) {{I(Gamma = 1) w1/h11 + I(Gamma = 2) w1/h11}
  wt1 <- w1h12*AItR
    
  ##  dtilNu(t) and Ytilu(t)
    
  dNtTilde <- dNt*{wt0 + wt1}

  return( list("dNtTilde" = dNtTilde,
               "Y" = Y,
               "wt0" = wt0, 
               "wt1" = wt1) )
}
