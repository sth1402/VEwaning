#' Estimation of Vaccine Efficacy Ove Time
#'
#' Implements a statistical framework for evaluating the efficacy of vaccines 
#'  based on a potential outcomes formulation.
#'
#' Note the infection time, U, can take 
#'   values NA or a value > L if the participant did not become infected.
#'  All other data must be complete.
#'  
#'
#' @param data A data.frame object containing all relevant data.
#'
#' @param L A numeric object. The analysis time.
#'
#' @param ... Ignored. Used only to require named inputs.
#'
#' @param modelGam1 A formula object. The coxph model for Gamma = 1.
#'   The LHS is set as the appropriate Surv() object internally. If a LHS
#'   is provided, it is ignored.
#'
#' @param modelGam2 A formula object. The coxph model for Gamma = 2.
#'   The LHS is set as the appropriate Surv() object internally. If a LHS
#'   is provided, it is ignored.
#'
#' @param modelEntry A formula object. The coxph model for entry times.
#'   The LHS is set as the appropriate Surv() object internally. If a LHS
#'   is provided, it is ignored.
#'
#' @param modelPsiGam1 A formula object. The logistic model for vaccination
#'   for participants with Gamma = 1.
#'   If a LHS is provided, it is ignored.
#'
#' @param modelPsiGam2 A formula object. The logistic model for vaccination
#'   for participants with Gamma = 2.
#'   If a LHS is provided, it is ignored.
#'
#' @param gFunc A character object. The model of infection rates. Must be one
#'   of {'lin', 'piece', 'splin', 'spcub'} for the linear,
#'   piecewise constant, linear spline, and cubic spline models respectively
#'
#' @param lag A numeric object. The lag time between the initial vaccine
#'   dose and full efficacy.
#'
#' @param v A numeric vector. The knots or cut-offs to be used by gFunc.
#'   If gFunc = 'lin', this input is ignored. For 'splin' and 'spcub', this
#'   is the knots of the spline on (0,L). For 'piece', v is the cut-offs on 
#'   (0,L). Note that this input should not include the extremes 0 and L.
#'
#' @param minWgt A numeric object. If not NULL, the minimum non-zero value a 
#'   weight can have, i.e., weight = max(minWgt, weight). If NULL, no
#'   truncation of weights is performed.
#'
#' @param maxWgt A numeric object. If not NULL, the maximum value a 
#'   weight can have, i.e., weight = min(maxWgt, weight). If NULL, no
#'   truncation of weights is performed.
#'
#' @param txName A character object. The header of the column of data 
#'   containing the treatment variable. Default value is 'A'.
#'   Treatment must be coded as 0/1, where 1 indicates that participant
#'   was vaccinated; 0 otherwise.
#'
#' @param infectionTime A character object. The header of the column of data 
#'   containing the time of infection on the scale of the calendar time. 
#'   Default value is 'U'.
#'
#' @param entryTime A character object. The header of the column of data 
#'   containing the time of entry into the study on the scale of the
#'   calendar time. Default value is 'E'.
#'
#' @param Gamma A character object. The header of the column of data 
#'   containing the category for the unblinding dynamic. 
#'   Default value is 'Gam'.
#'   Data must be 0/1/2, where 0 indicates infection occurs before requested/
#'   offered unblinding; 1 indicates unblinding was requested by participant
#'   prior to the commencement of participant decision clinic visits; 
#'   and 2 indicates that unblinding occurred after the commencement of
#'   participant decision clinic visits
#'
#' @param unblindTime A character object. The header of the column of data 
#'   containing the time to requested unblinding, participant  decision
#'   clinic visit/requested unblinding, or infection, whichever comes first.
#'   Default value is 'R'.
#'
#' @param vaccinated A character object. The header of the column of data 
#'   containing the indicator of vaccination, where 1 if participant is
#'   vaccinated; 0 otherwise.
#'   Default value is 'Psi'.
#'
#' @returns A list object
#'   \item{theta}{A vector object containing the estimated theta parameters.}
#'   \item{cov}{The covariance estimated using the sandwich estimator.}
#'   \item{SE}{The standard error estimated using the sandwich estimator.}
#'
#' @include unblindFit.R entryFit.R psiFit.R
#'
#' @references Tsiatis, A. A. and Davidian, M. (2021) Estimating Vaccine
#'   Efficacy Over Time After a Randomized Study is Unblinded. Submitted.
#'
#' @export
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
#' @import survival
#' @importFrom stats as.formula coef predict update.formula
veWaning <- function(data, 
                     L, 
                     ..., 
                     lag = 0.0,
                     modelGam1 = NULL,
                     modelGam2 = NULL,
                     modelEntry = NULL,
                     modelPsiGam1 = NULL,
                     modelPsiGam2 = NULL,
                     gFunc = NULL,
                     v = NULL,
                     minWgt = NULL,
                     maxWgt = NULL,
                     txName = "A",
                     infectionTime = "U", 
                     entryTime = "E", 
                     Gamma = "Gam", 
                     unblindTime = "R",
                     vaccinated = "Psi") {

  # list to store the column names of relevant data
  dObj <- list("Gam" = Gamma, 
               "R" = unblindTime, 
               "E" = entryTime,
               "U" = infectionTime,
               "A" = txName,
               "Psi" = vaccinated)

  # ensure that the names exist in the provided data.frame
  tryCatch(expr = data[,c(dObj$Gam, dObj$R, dObj$E, dObj$U, dObj$A, dObj$Psi)],
           error = function(e) {
                     message("unable to identify key components of input data")
                     stop(e$message, call. = FALSE)
                   })

  # ensure that the Gamma is an integer
  iGam <- as.integer(x = round(x = data[,dObj$Gam], digits = 0L))

  # ensure that Gamma is one of (0,1,2)
  if (any(!{iGam %in% c(0L,1L,2L)})) {
    stop("unrecognized Gamma values", call. = FALSE)
  }
  data[,dObj$Gam] <- iGam

  # ensure that treatment is integer
  iA <- as.integer(x = round(x = data[,dObj$A], digits = 0L))

  # ensure that treatment is one of (0,1)
  if (any(!{iA %in% c(0L,1L)})) {
    stop("unrecognized treatment values", call. = FALSE)
  }
  data[,dObj$A] <- iA

  # ensure that Psi is integer
  iPsi <- as.integer(x = round(x = data[,dObj$Psi], digits = 0L))

  # ensure that Psi is one of (0,1)
  if (any(!{iPsi %in% c(0L,1L)})) {
    stop("unrecognized indicator of vaccination values", call. = FALSE)
  }
  data[,dObj$Psi] <- iPsi

  # ensure gFunc is appropriately specified
  if (!is.null(x = gFunc)) {

    # convert to lower
    gFunc <- tolower(x = gFunc)

    # must be one of (lin, piece, splin, spcub)
    if (!{gFunc %in% c('lin', 'piece', 'splin', 'spcub')}) {
      stop("gFunc is not recognized", call. = FALSE)
    }

    # if it is (piece, splin, spcub) knots or cutoff must be provided
    if (gFunc %in% c('piece', 'splin', 'spcub')) {

      if (is.null(x = v)) {
        stop("v must be provided for the selected gFunc", call. = FALSE)
      }

      if (gFunc %in% 'piece') {
        # piecewise needs lower and upper bounds to be 0, L respectively
        if (min(v) > 0.0) v <- c(0.0, v)
        if (max(v) < L) v <- c(v,L)
        nTheta <- length(x = v) - 1L
      } else if (gFunc %in% 'splin') {
        # linear spline has two additional thetas
        nTheta <- length(x = v) + 2L
      } else if (gFunc %in% 'spcub') {
        # cubic spline has four additional thetas
        nTheta <- length(x = v) + 4L
      }
    } else {
      nTheta <- 2L
    }
  } else {
    gFunc <- 'lin'
    nTheta <- 2L
  }

  # testing to ensure that none or all models are provided
  tst <- c(is.null(x = modelGam1),
           is.null(x = modelGam2), is.null(x = modelEntry),
           is.null(x = modelPsiGam1), is.null(x = modelPsiGam2))

  if (!all(tst) && !all(!tst)) {
    stop("insufficient models provided", call. = FALSE)
  }

  # determine indicator of infection by time L
  data[is.na(x = data[,dObj$U]),dObj$U] <- L + 10.0
  delta <- as.integer(x = data[,dObj$U] < L)

  # limit vaccination times to L
  data[,dObj$U] <- data[,dObj$U] * {delta == 1L} + 
                   L * {delta == 0L}

  # times of infection for those not unblinded
  timesB <- sort(x = unique(x = data[delta == 1L & data[,dObj$Gam] == 0L, dObj$U]))
  nTimesB <- length(x = timesB)

  # times of infection for unblinded
  timesU <- sort(x = unique(x = data[delta == 1L & data[,dObj$Gam] != 0L, dObj$U]))
  nTimesU <- length(x = timesU)

  # earliest time that unblinding was requested
  TP <- min(data[data[,dObj$Gam] != 0L, dObj$R]) 

  # latest time that unblinding occured
  TC <- max(data[data[,dObj$Gam] != 0L, dObj$R]) 

  # earliest time that unblinding was scheduled
  TU <- min(data[data[,dObj$Gam] == 2L, dObj$R])

  if (any(!tst)) {

    KRfR <- unblindFit(data = data, 
                       times = timesB,  
                       TP = TP,  
                       TU = TU,
                       TC = TC,  
                       modelGam1 = modelGam1,  
                       modelGam2 = modelGam2,  
                       dObj = dObj)

    fEX <- entryFit(data = data, modelEntry = modelEntry, dObj = dObj)

    Psiprob <- psiFit(data = data, 
                      modelPsi1 = modelPsiGam1, 
                      modelPsi2 = modelPsiGam2, 
                      dObj = dObj)

    weightOne <- FALSE

  } else {

    n <- nrow(x = data)

    KRfR <- list("KR.0.stab" = matrix(data = 1.0, nrow = n, ncol = nTimesB),
                 "KR.1.stab" = matrix(data = 1.0, nrow = n, ncol = nTimesB),
                 "fR1.0.stab" = rep(x = 1.0, times = n),
                 "fR1.1.stab" = rep(x = 1.0, times = n),
                 "fR2.0.stab" = rep(x = 1.0, times = n),
                 "fR2.1.stab" = rep(x = 1.0, times = n))

    fEX <- rep(x = 1.0, times = n)

    Psiprob <- list("pPsi1.stab" = rep(x = 1.0, times = n),
                    "pPsi2.stab" = rep(x = 1.0, times = n))

    weightOne <- TRUE
  }

  if (is.null(x = minWgt)) minWgt <- 0.0
  if (is.null(x = maxWgt)) maxWgt <- Inf

  out <- esttheta(data = data,
                  timesB = timesB,
                  timesU = timesU,
                  KRfR = KRfR,
                  fEX = fEX,
                  Psiprob = Psiprob,
                  lag = lag,
                  gFunc = gFunc,
                  v = v,
                  nTheta = nTheta,
                  dObj = dObj,
                  delta = delta,
                  minWgt = minWgt,
                  maxWgt = maxWgt)


  allTau <- {data[,dObj$A] == 1L}*{L - data[,dObj$E]} +
            {data[,dObj$A] == 0L}*{data[,dObj$Psi] == 1L}*{L - data[,dObj$R]}

  maxTau <- max(allTau)
  if (maxTau < lag) {
    message("no participant has been vaccinated longer than the lag time")
    maxTau <- NULL
  }

  attr(out, "gFunc") <- gFunc
  attr(out, "maxTau") <- maxTau
  attr(out, "lag") <- lag
  attr(out, 'v') <- v

  class(out) <- "VEwaning"

  return( out )
    
}


