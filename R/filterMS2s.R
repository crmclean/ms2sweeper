#' @title filterMS2s
#'
#' @description This function is designed to apply user input filters to filter
#' ms2 peaks of the data. Currently this funtion filters by scan count and by
#' minimum intensity value.
#'
#' @param sweeperObj - Sweeper Object containing MS2 spectra.
#' @param minScanCount - Minimum number of scan containing a peak.
#' @param minIntensity - Minimum intensity value required to retain a peak.
#' @param minScanScore - Minimum possible scan score to retain.
#'
#' @return returns the sweeperObj with the filttered ms2 spectra by the
#' user supplied heuristic values.
#' @export
filterMS2s <- function(sweeperObj, minScanCount = 0, minIntensity = 10,
                       minScanScore = 10) {

    ms2Scans <- getMs2Pure(sweeperObj)
    if(length(ms2Scans) == 0) {
        stop("Pure MS2 scans must be loaded to sweeper object in order to filter them.")
    }


    for(i in seq_along(ms2Scans)) {
        curScan <- ms2Scans[[i]]
        curScan <- curScan[curScan$intMean > minIntensity,]

        curScan <- curScan[curScan$peakScore > minScanScore,]


        if(max(curScan$size) > 1) {
            curScan <- curScan[curScan$size > minScanCount,]
        }
        ms2Scans[[i]] <- curScan
    }

    sweeperObj <- setMs2Pure(sweeperObj, ms2Scans)
    return(sweeperObj)

}
