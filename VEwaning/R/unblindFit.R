pzero <- function(object, newdata) {
  temp <- predict(object = object, 
                  newdata = newdata, 
                  type = "lp", 
                  reference = "sample") + 
          sum(coef(object = object) * object$means, na.rm=TRUE)
  return( exp(x = temp) )
}


######################################################################

##  Fit Cox models for unblinding hazards for Gam = 1 and 2 and
##  calculate survival probabilities K_R at U's for all individuals
##  with Delta=1 whose U's fall in [TP,TU), Gam=1 and [TU,TC), Gam=2.
##  Also calculate the survival probabilities at the mean X values to
##  form stabilized weights and get the f_R1 and f_R2 at the R values

##  merge() used below will order the KR1 and KR2 data frames by ID
##  number in original data set
#
#' @importFrom stats model.response
#' @importFrom stats model.frame
#' @import survival
#'
# @param data A data.frame object. All data required for the analysis.
#
# @param infectionTimes A vector object. Infection times.
#
# @param TP A numeric object. Time at which Pfizer granted EUA.
#
# @param TU A numeric object. Time at which Participant Decision Clinic Visits
#   commenced.
#
# @param TC A numeric object. Time at which Participant Decision Clinic Visits
#   ended.
#
# @param modelGam1 A formula object. The coxph model for Gamma = 1. Note
#   that the LHS is taken as a Surv() object as defined by package survival.
#
# @param modelGam2 A formula object. The coxph model for Gamma = 2. Note
#   that the LHS is taken as a Surv() object as defined by package survival.
#
# @param dObj A list object. The column headers of data containing required
#   variables
#
unblindFit <- function(data, 
                       times,  
                       TP, 
                       TU, 
                       TC,  
                       modelGam1,  
                       modelGam2,  
                       dObj) {

  n <- nrow(x = data)

  # remove LHS of each model
  modelGam1 <- update.formula(old = modelGam1, new = NULL ~ .)
  modelGam2 <- update.formula(old = modelGam2, new = NULL ~ .)

  # ensure that model frame can be constructed from provided formula
  mFrame1 <- tryCatch(expr = stats::model.frame(formula = modelGam1, 
                                                data = data),
                      error = function(e) {
                                message("unable to obtain model.frame")
                                stop(e$message, call. = FALSE)
                              })

  mFrame2 <- tryCatch(expr = stats::model.frame(formula = modelGam2, 
                                                data = data),
                      error = function(e) {
                                message("unable to obtain model.frame")
                                stop(e$message, call. = FALSE)
                              })

  # keep only the covariates used in the two models
  covs <- sort(x = unique(x = c(colnames(x = mFrame1), colnames(x = mFrame2))))
  nCovs <- length(x = covs)

  # add Gamma and R data
  data <- cbind(data[,covs], data[,dObj$Gam], data[,dObj$R])
  colnames(x = data) <- c(covs, dObj$Gam, dObj$R)
  timeCol <- ncol(x = data)

  # update formulae to include Surv() LHS
  survForm <- paste0('Surv(',dObj$R,',1*{',dObj$Gam,' == 1L})')
  modelGam1 <- paste0(survForm, "~", as.character(x = modelGam1)[2L])
  modelGam1 <- as.formula(object = modelGam1)

  survForm <- paste0('Surv(',dObj$R,',1*{',dObj$Gam,' == 2L})')
  modelGam2 <- paste0(survForm, "~", as.character(x = modelGam2)[2L])
  modelGam2 <- as.formula(object = modelGam2)

  #  Fit Cox models for R in [TP,TU), Gam=1 and [TU,TC), Gam=2
    
  gam1Fit <- tryCatch(expr = survival::coxph(formula = modelGam1, data = data),
                      error = function(e) {
                                message("unable to Cox model for R in [TP,TU) Gam = 1")
                                stop(e$message, call. = FALSE)
                              })

  if (any(is.na(x = gam1Fit$coef))) {
    stop("fit of Cox model for R in [TP,TU) Gam = 1 results in NA coefficients",
         call. = FALSE)
  }

  gam2Fit <- tryCatch(expr = survival::coxph(formula = modelGam2, data = data),
                      error = function(e) {
                                message("unable to fit Cox model for R in [TP,TU) Gam = 2")
                                stop(e$message, call. = FALSE)
                              })

  if (any(is.na(x = gam2Fit$coef))) {
    stop("fit of Cox model for R in [TP,TU) Gam = 2 results in NA coefficients",
         call. = FALSE)
  }

  # means of all covariates. This includes Gamma and R, but they will be
  # reset prior to use
  meansAll <- colMeans(x = data)
  means0 <- colMeans(x = data[data[,dObj$A] == 0L,,drop = FALSE])
  means1 <- colMeans(x = data[data[,dObj$A] == 1L,,drop = FALSE])

  ##  Survival probabilities for Gam=1 at U's in [TP,TU) -- last row of
  ##  KR1.all will be at mean X, A values, so extract separately

  gam1 <- {times >= {TP-1e-8}} & {times < TU}

  # {nTimesGam1}
  timesGam1 <- data.frame(times[gam1])
  colnames(x = timesGam1) <- dObj$R

  # {n+2 x nM}
  newdata <- rbind(data, means0, means1)
  newdata[,dObj$Gam] <- 1L

  newdata <- merge(x = newdata[,-timeCol], y = timesGam1)

  # {n+2 x nTimesGam1}
  predsurvR1 <- matrix(data = exp(x = -predict(object = gam1Fit,
                                               newdata = newdata,
                                               type = "expected")),
                       nrow = n + 2L, 
                       ncol = nrow(x = timesGam1))

  # {nTimesGam1 x n+2}
  predsurvR1 <- t(x = predsurvR1)

  #  Get survival probabilities at TU for all individuals 

  # {n+2 x nM}
  newdata <- rbind(data, means0, means1)
  newdata[,dObj$Gam] <- 1L
  newdata[,dObj$R] <- TU

  # {n+2}
  predsurvTU <- drop(x = exp(x = -predict(object = gam1Fit,
                                          newdata = newdata,
                                          type = "expected")))
        
  #  Survival probabilities for Gam=2 at U's in [TU,TC) -- last row of
  #  KR2.all will be at mean X values, extract separately
  
  gam2 <- times >= {TU-1e-8} & times < TC

  # {nTimesGam2}
  timesGam2 <- data.frame(times[gam2])
  colnames(x = timesGam2) <- dObj$R

  newdata <- rbind(data, meansAll, meansAll)
  newdata[,dObj$Gam] <- 2L

  newdata <- merge(x = newdata[,-timeCol], y = timesGam2)

  # {n+2 x nTimesGam2}
  predsurvR2 <- matrix(data = exp(x = -predict(object = gam2Fit,
                                               newdata = newdata,
                                               type = "expected")),
                       nrow = n + 2L, ncol = nrow(x = timesGam2))

  ##  Multiply these by the prob at TU
  # {nTimesGam2 x n+2}
  KR12 <- t(x = predsurvTU*predsurvR2)

  ##  Form stabilized quantities for stabilized weights

  ##  A=0 

  # {n x nTimesGam1}
  KR1.0.stab <- 1.0 / t(x = predsurvR1[,1L:n,drop = FALSE] / predsurvR1[,n + 1L])

  # {n x nTimesGam2}
  KR12.0.stab <- 1.0 / t(x = KR12[,1L:n,drop = FALSE] / KR12[,n + 1L])

  ##  A=1

  # {n x nTimesGam1}
  KR1.1.stab <- 1.0 / t(x = predsurvR1[,1L:n,drop = FALSE] / predsurvR1[,n + 2L])

  # {n x nTimesGam2}
  KR12.1.stab <- 1.0 / t(x = KR12[,1L:n,drop = FALSE] / KR12[,n + 2L])
    
  ##  Form matrix of stabilized KR for U's <=TC   

  gam0 <- times < TP
  timesGam0 <- times[gam0]

  # {n x nTimesGam0}
  KR0 <- matrix(data = 1.0, nrow = n, ncol = length(x = timesGam0))

  ##  These are n x lUd0+lUd1+lUd2 matrices with columns for each U<TC
  # {n x nTimesGam0+nTimesGam1+nTimesGam2}

  KR.0.stab <- cbind(KR0, KR1.0.stab, KR12.0.stab)
  KR.1.stab <- cbind(KR0, KR1.1.stab, KR12.1.stab)

  ##  Form stabilized f_R densities evaluated at the R values for Gam =
  ##  0, 1, and 2, respectively -- first get KR -- because predicted
  ##  KR1r=1 for r < TP, KR2r=1 for r<TU, and KR1r=approx surv prob at
  ##  TU for r>TU, all we have to do to get KR matrix with correct
  ##  survival probs for each range is multiply

  newdata <- data

  newdata[,dObj$A] <- 0L
  KR1r.0 <- exp(x = -predict(object = gam1Fit,
                             newdata = newdata,
                             type = "expected"))
  eR1.0 <- pzero(object = gam1Fit, newdata = newdata)

  newdata[,dObj$A] <- 1L
  KR1r.1 <- exp(x = -predict(object = gam1Fit,
                             newdata = newdata,
                             type = "expected"))
  eR1.1 <- pzero(object = gam1Fit, newdata = newdata)


  KR2r <- exp(x = -predict(object = gam2Fit,
                           newdata = data,
                           type = "expected"))
  eR2 <- pzero(object = gam2Fit, newdata = data)

  KRr.0 <- KR1r.0*KR2r
  KRr.1 <- KR1r.1*KR2r

  ##  For stabilized weights must get the same thing for the means of X
  ##  separately by A

  newdata <- data
  for (i in 1L:nCovs) {
    newdata[,i] <- means0[i]
  }

  KR1r.mean.0 <- exp(x = -predict(object = gam1Fit,
                                  newdata = newdata,
                                  type = "expected"))
  eR1.mean.0 <- pzero(object = gam1Fit, newdata = as.data.frame(t(x = means0)))

  for (i in 1L:nCovs) {
    newdata[,i] <- means1[i]
  }

  KR1r.mean.1 <- exp(x = -predict(object = gam1Fit,
                                  newdata = newdata,
                                  type = "expected"))
  eR1.mean.1 <- pzero(object = gam1Fit, newdata = as.data.frame(t(x = means1)))

  for (i in 1L:nCovs) {
    newdata[,i] <- meansAll[i]
  }

  KR2r.mean <- exp(x = -predict(object = gam2Fit,
                                newdata = newdata,
                                type = "expected"))
  eR2.mean <- pzero(object = gam2Fit, newdata = as.data.frame(t(x = meansAll)))


  KRr.mean.0 <- KR1r.mean.0*KR2r.mean
  KRr.mean.1 <- KR1r.mean.1*KR2r.mean

  ##  Stabilized densities at each R for each A, n x 1 vectors

  fR1.0.stab <- {eR1.mean.0*KRr.mean.0} / {eR1.0*KRr.0}
  fR1.1.stab <- {eR1.mean.1*KRr.mean.1} / {eR1.1*KRr.1}
  fR2.0.stab <- {eR2.mean*KRr.mean.0} / {eR2*KRr.0}
  fR2.1.stab <- {eR2.mean*KRr.mean.1} / {eR2*KRr.1}
    
  return( list("KR.0.stab" = KR.0.stab,
               "KR.1.stab" = KR.1.stab,
               "fR1.0.stab" = fR1.0.stab,
               "fR1.1.stab" = fR1.1.stab,
               "fR2.0.stab" = fR2.0.stab,
               "fR2.1.stab" = fR2.1.stab) )
}
