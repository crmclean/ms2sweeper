#' @title interpolateIntensity
#'
#' @description This funciton is meant to interpolate the proper intensity value for the MS2
#' between the two adjacent peaks. This is done in attempt to correct the spectra from convolution
#' due to noise of features.
interpolateIntensity <- function(ms2Matches,
                                 purityTable,
                                 curMs2) {

  featureRt <- ms2Matches[[curMs2]]$ms2rt[1]


  y = purityTable$intensity
  x = purityTable$scanTime

  lineCoefs <- coef(lm(y ~ x))

  return(lineCoefs)

}
