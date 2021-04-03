##  Fit a Cox model to the entry times as a function of X and return a
##  n x 1 vector of stabilized densities evaluated at the observed E, X

# @param data A data.frame object. All data required for the analysis.
#
# @param modelEntry A formula object. The coxph model for Entry time. Note
#   that the LHS is taken to be a Surv() object as defined by package survival.
#
# @param dObj A list object. The column headers of data containing required
#   variables
#
#' @importFrom stats model.frame
#' @importFrom stats model.response
#' @import survival
entryFit <- function(data, modelEntry, dObj) {

  # create internal naming convention for Status
  while (TRUE) {
    statusName = paste(sample(x = letters, size = 10, replace = TRUE), 
                       collapse = "")
    if (statusName %in% colnames(x = data)) next
    break
  }

  # remove LHS of model
  modelEntry <- update.formula(old = modelEntry, new = NULL ~ .)

  mFrame1 <- tryCatch(expr = stats::model.frame(formula = modelEntry, 
                                                data = data),
                      error = function(e) {
                                message("unable to obtain model.frame")
                                stop(e$message, call. = FALSE)
                              })
    
  # keep only the covariates used in the model
  covs <- colnames(x = mFrame1)
  nCovs <- length(x = covs)

  # redefine the data and formula to correspond to internal naming
  # convention
  data <- cbind(data[,covs], 1L, data[,dObj$E])
  colnames(x = data) <- c(covs, statusName, dObj$E)
  timeCol <- ncol(x = data)

  # update formula with Surv() object in LHS
  entryForm <- paste0('Surv(',dObj$E,',',statusName,')')
  modelEntry <- paste0(entryForm, "~", as.character(x = modelEntry)[2L])
  modelEntry <- as.formula(object = modelEntry)

  EObj <- tryCatch(expr = survival::coxph(formula = modelEntry, data = data),
                   error = function(e) {
                             message("unable to fit Entry time cox model")
                             stop(e$message, call. = FALSE)
                           })

  if (any(is.na(x = EObj$coef))) {
    stop("fit of Cox model for Entry time results in NA coefficients",
         call. = FALSE)
  }

  ##  Predicted survival probabilities at each E

  SE <- exp(x = -predict(object = EObj, newdata = data, type = "expected"))
  eE <- pzero(object = EObj, newdata = data)

  ##  Predicted survival probabilities at each E for mean of X    

  newdata <- data
  means <- colMeans(x = data)
  for (i in 1L:nCovs) {
    newdata[,i] <- means[i]
  }
  
  SE.mean <- exp(x = -predict(object = EObj, 
                              newdata = newdata, 
                              type = "expected"))

  eE.mean <- pzero(object = EObj, newdata = newdata[1L,])

  ##  Stabilized density for each individual

  fE.stab <- {eE.mean*SE.mean} / {eE*SE}
    
  return( fE.stab )
}
