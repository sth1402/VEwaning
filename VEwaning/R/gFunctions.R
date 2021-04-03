gFunction <- function(gFunc, u, theta, knts) {

  if (gFunc == 'piece') {
    return( gPiecewise(u = u, theta = theta, v = knts) )
  } else if (gFunc == 'lin') {
    return( gLinear(u = u, theta = theta) )
  } else if (gFunc == 'splin') {
    return( gSplineLinear(u = u, theta = theta, knots = knts) )
  } else if (gFunc == 'spcub') {
    return( gSplineCubic(u = u, theta = theta, knots = knts) )
  } else {
    stop("unrecognized gFunc", call. = FALSE)
  }
}

gPiecewise <- function(u, theta, ..., v) {

  if (!is.matrix(x = u)) u <- matrix(data = u, ncol = 1L)

  gu <- matrix(data = 0.0, nrow = nrow(x = u), ncol = ncol(x = u))
  gutheta <- list()

  gutheta[[ 1L ]] <- matrix(data = 1.0, nrow = nrow(x = u), ncol = ncol(x = u))

  i <- 3L
  while (i <= length(x = v)) {
    gutheta[[ i - 1L ]] <- 1.0*{{u > v[i-1L]} & {u <= v[i]}}
    gu <- gu + theta[i-1L]*gutheta[[ i - 1L ]]
    i <- i + 1L
  }

  return( list("gu" = gu, "gutheta" = gutheta) )

}   

gLinear <- function(u, theta, ...) {

  if (!is.matrix(x = u)) u <- matrix(data = u, ncol = 1L)

  gu <- matrix(data = 0.0, nrow = nrow(x = u), ncol = ncol(x = u))
  gutheta <- list()

  gutheta[[ 1L ]] <- matrix(data = 1.0, nrow = nrow(x = u), ncol = ncol(x = u))

  gu <- gu + theta[2L]*u
  gutheta[[ 2L ]] <- u

  return( list("gu" = gu, "gutheta" = gutheta) )

}   

gSplineLinear <- function(u, theta, ..., knots) {

  if (!is.matrix(x = u)) u <- matrix(data = u, ncol = 1L)

  gu <- matrix(data = 0.0, nrow = nrow(x = u), ncol = ncol(x = u))
  gutheta <- list()

  gutheta[[ 1L ]] <- matrix(data = 1.0, nrow = nrow(x = u), ncol = ncol(x = u))

  gu <- theta[2L]*u
  gutheta[[ 2L ]] <- u

  for (k in 1L:length(x = knots)) {
    temp <- {u-knots[k]}
    temp[u <= knots[k]] <- 0.0
    gu <- gu + theta[k+2L]*temp
    gutheta[[ k+2L ]] <- temp
  }

  return( list("gu" = gu, "gutheta" = gutheta) )

}

gSplineCubic <- function(u, theta, ..., knots) {

  if (!is.matrix(x = u)) u <- matrix(data = u, ncol = 1L)

  gu <- matrix(data = 0.0, nrow = nrow(x = u), ncol = ncol(x = u))
  gutheta <- list()

  gutheta[[ 1L ]] <- matrix(data = 1.0, nrow = nrow(x = u), ncol = ncol(x = u))

  gutheta[[ 2L ]] <- u
  gu <- gu + theta[2L]*gutheta[[ 2L ]]

  gutheta[[ 3L ]] <- u^2L
  gu <- gu + theta[3L]*gutheta[[ 3L ]]

  gutheta[[ 4L ]] <- u^3L
  gu <- gu + theta[4L]*gutheta[[ 4L ]]

  for (k in 1L:length(x = knots)) {
    temp <- {u-knots[k]}^3L*{u > knots[k]}
    gu <- gu + theta[k+4L]*temp
    gutheta[[ k+4L ]] <- temp
  }

  return( list("gu" = gu, "gutheta" = gutheta) )

}
