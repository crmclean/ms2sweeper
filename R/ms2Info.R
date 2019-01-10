#' @title ms2Info
#'
#' @description This function is used to subset entered peak table for all fragments into only values
#' related to MS2 observations.
#'
#' @param MSn_dt - list of tables containing all fragments observed
ms2Info <- function(MSn_dt) {
  ms2Level <- MSn_dt$msLevel == 2
  ms2mz <- MSn_dt$mz[ms2Level]
  ms2rt <- MSn_dt$rt[ms2Level]
  ms2intensity <- MSn_dt$intensity[ms2Level]
  ms2parentPeaks <- MSn_dt$MSnParentPeakID[ms2Level]
  ms2sample <- MSn_dt$dataFile[ms2Level]
  return(list(ms2mz = ms2mz,
              ms2rt = ms2rt,
              ms2intensity = ms2intensity,
              ms2parentPeaks = ms2parentPeaks,
              ms2sample = ms2sample))
}
