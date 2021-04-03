##  Fit logistic regression model for Psi with A=0 for Gam = 1, 2
##  and get stabilized estimated probabilities for each individual's Gam

#' @importFrom stats glm
#' @importFrom stats predict.glm
psiFit <- function(data, modelPsi1, modelPsi2, dObj) {

  dat0 <- data[ {data[,dObj$A] == 0L} & {data[,dObj$Gam] != 0L},]

  grp1 <- {data[,dObj$A] == 0L} & {data[,dObj$Gam] == 1L}
  means1 <- colMeans(x = data[grp1,])

  grp2 <- {data[,dObj$A] == 0L} & {data[,dObj$Gam] == 2L}
  means2 <- colMeans(x = data[grp2,])

  # remove LHS of formulae
  modelPsi1 <- update.formula(old = modelPsi1, new = NULL ~ .)
  modelPsi2 <- update.formula(old = modelPsi2, new = NULL ~ .)

  # update formulae to include Psi
  modelPsi1 <- paste0(dObj$Psi, "~", as.character(x = modelPsi1)[2L])
  modelPsi1 <- as.formula(object = modelPsi1)

  modelPsi2 <- paste0(dObj$Psi, "~", as.character(x = modelPsi2)[2L])
  modelPsi2 <- as.formula(object = modelPsi2)

  ##  Fit logistic regression model only to these data
    
  logist1 <- tryCatch(expr = glm(formula = modelPsi1,
                                 data = data[grp1,],
                                 family = 'binomial'),
                      error = function(e) {
                                message("unable to obtain fit for psi1")
                                stop(e$message, call. = FALSE)
                              })
  if (any(is.na(x = logist1$coef))) {
    stop("glm fit returns NA coefficients for modelPsiGam1",
         call. = FALSE)
  }

  logist2 <- tryCatch(expr = glm(formula = modelPsi2,
                                 data = data[grp2,],
                                 family = 'binomial'),
                      error = function(e) {
                                message("unable to obtain fit for psi2")
                                stop(e$message, call. = FALSE)
                              })

  if (any(is.na(x = logist2$coef))) {
    stop("glm fit returns NA coefficients for modePsiGam2",
         call. = FALSE)
  }
    
  ##  Get predicted probabilities for Gam = 1, 2 for everyone     
  pPsi1 <- stats::predict.glm(object = logist1, 
                              newdata = data, 
                              type = "response")

  pPsi2 <- stats::predict.glm(object = logist2, 
                              newdata = data, 
                              type = "response")

  ##  Get predicted probabilities at the means for Gam = 1 and Gam = 2
  pPsi1.mean <- stats::predict.glm(object = logist1, 
                                   newdata = as.data.frame(t(x = means1)), 
                                   type = "response")

  pPsi2.mean <- stats::predict.glm(object = logist2, 
                                   newdata = as.data.frame(t(x = means2)), 
                                   type = "response")

  ##  Get stabilized probabilities
    
  pPsi1.stab <- pPsi1.mean / pPsi1
  pPsi2.stab <- pPsi2.mean / pPsi2

  return( list("pPsi1.stab" = pPsi1.stab,
               "pPsi2.stab" = pPsi2.stab) )
   
}
