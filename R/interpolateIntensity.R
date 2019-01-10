#' @title interpolateIntensity
#'
#' @description This funciton is meant to interpolate the proper intensity value for the MS2
#' between the two adjacent peaks. This is done in attempt to correct the spectra from convolution
#' due to noise of features.
#'
#' @param purityTable - Table of parent ion features observed over time.
#'
#' @return coefficients from running a regression on the parent ions over
#' the time window where the peak is formed.
interpolateIntensity <- function(purityTable) {

  y = purityTable$intensity
  x = purityTable$scanTime

  lineCoefs <- stats::coef(stats::lm(y ~ x))

  return(lineCoefs)

}
