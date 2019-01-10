#' @title check13c
#'
#' @param windowData - the mz intensity table making up the observations of a
#' given mass spectra.
#' @param ppm - error window used to identify 13C peaks
#'
#' @description This function checks if the scan with the parent ion has a
#' 13C feature. If it does, it remove the 13C feature from the scan window and
#' combines the purity. The new table tells us if we need to fix the ms2 should
#' be corrected for contamination.
#'
#' @return The windowData table with an additional column used to indicate
#' if the data has a 13C peak that could be contamination the ms2 scan.
check13c <- function(windowData, ppm) {

    windowData$has13C <- F
    if(nrow(windowData) == 1) {
        return(windowData)
    }

    if(abs((windowData$mzs[1] + 1) - windowData$mzs[2])/windowData$mzs[1] * 10^6 < ppm*2) {

        newPurityScore <- sum(windowData$purity[c(1,2)])
        windowData$purity[c(1,2)] <- newPurityScore
        windowData$intensity[c(1,2)] <- sum(windowData$intensity[c(1,2)])
        windowData <- windowData[-2,]
        windowData$has13C[1] <- T

    }

    return(windowData)
}
