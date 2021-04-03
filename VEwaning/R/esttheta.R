#' @include pre_dNYb.R pre_dNYu.R
esttheta <- function(data, 
                     timesB, 
                     timesU,
                     KRfR, 
                     fEX,  
                     Psiprob,  
                     lag, 
                     gFunc,
                     v,  
                     nTheta,
                     dObj,
                     delta,
                     minWgt,
                     maxWgt) {

  ##  Initialize -- starting value 

  theta <- rep(x = 0.0, times = nTheta)
  score <- rep(x = 1.0, times = nTheta)

  bmax <- Inf
  tol <- 1e-5
  imax <- 20L
  iter <- 1L

  # components of estimators that do not depend on theta
  preb <- pre_dNYb(data= data, 
                   times = timesB, 
                   KRfR = KRfR, 
                   fEX = fEX, 
                   lag = lag, 
                   dObj = dObj,
                   delta = delta,
                   minWgt = minWgt,
                   maxWgt = maxWgt)

  preu <- pre_dNYu(data= data, 
                   times = timesU, 
                   KRfR = KRfR, 
                   fEX = fEX, 
                   Psiprob = Psiprob,
                   lag = lag, 
                   dObj = dObj,
                   delta = delta,
                   minWgt = minWgt,
                   maxWgt = maxWgt)


  while (iter < imax && bmax > tol) {
  
    ## Get blinded contribution to estimating equation, gradient, 
    ## influence function
        
    eb <- estimateb(data = data,  
                    times = timesB,  
                    preb = preb,  
                    lag = lag, 
                    gFunc = gFunc,  
                    theta = theta,  
                    v = v,
                    dObj = dObj)

    ## Get blinded contribution to estimating equation, gradient, 
    ## influence function

    eu <- estimateu(data = data,  
                    times = timesU,  
                    preu = preu,  
                    lag = lag, 
                    gFunc = gFunc,  
                    theta = theta,  
                    v = v,
                    dObj = dObj)

    ##  Update theta        

    score <- eb$estb + eu$estu

    theta <- theta + solve(a = eb$gradb+eu$gradu, b = score)

    iter <- iter + 1L
    bmax <- max(abs(x = score))
  }

  ##  For final value of theta, get covariance matrix and compute estimated 
  ##  SEs -- we only use sandwich

  eb <- estimateb(data = data,  
                  times = timesB,  
                  preb = preb,  
                  lag = lag, 
                  gFunc = gFunc,  
                  theta = theta,  
                  v = v,
                  dObj = dObj)

  eu <- estimateu(data = data,  
                  times = timesU,  
                  preu = preu,  
                  lag = lag, 
                  gFunc = gFunc,  
                  theta = theta,  
                  v = v,
                  dObj = dObj)

  Covmodel <- solve(a = eb$gradb + eu$gradu)
  Covsand <- Covmodel %*% {eb$meatb + eu$meatu} %*% Covmodel
  SEsand <- sqrt(x = diag(x = Covsand))

  nms <- paste0('theta', 0L:{nTheta - 1L})

  names(x = theta) <- nms

  colnames(x = Covsand) <- nms
  rownames(x = Covsand) <- nms
  names(x = SEsand) <- nms


  return( list("theta" = theta,
               "cov" = Covsand,
               "SE" = SEsand) )
}
