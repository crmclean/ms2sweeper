#' @title parentPurity
#'
#' @description This is a function used to determine the purity of each parent
#' ion in making the ms2 spectra. The amuWindow defines the range of masses
#' that can enter the collision cell. Ideally, this will be used to help correct
#' patterns and deconvolute the ms2 spectra.
#'
#' @param sampleData - current sample of MS2 data being looked at.
#' @param curMz - m/z being queried within data.
#' @param parentMatches - index of sampleData that contain matches to curMz
#' @param allMasses - all mz values from sample raw data
#' @param allIntensities - all intensity values from sample raw data
#' @param scanTimes - all scans (rts) from sample raw data
#' @param ppm - mass error
#' @param isoWindow - a measure in dalton of the isolation error from the instrument
#' @param scanInterval - dynamic exclusion isolation window.
#'
#' @return List of parent m/z matches over time for each feature along with their purity.
parentPurity <- function(sampleData,
                         curMz,
                         parentMatches,
                         allMasses,
                         allIntensities,
                         scanTimes,
                         ppm,
                         isoWindow,
                         scanInterval) {

    parentRts <- sampleData$rt[parentMatches]
    parentIntensities <- sampleData$intensity[parentMatches]

    # Cleaning up each MS2 table ----------------------------------------------
    prevFile <- ""
    scansMassSpectra <- list()
    for(curMs2 in seq_along(parentRts)) {

        parentScan <- scanTimes == parentRts[curMs2]

        ## if there isn't an exact match, expand the rt window
        if(!any(parentScan)) {
            parentScan <- scanTimes > (parentRts[curMs2] - 0.5) &
                scanTimes < (parentRts[curMs2] + 0.5)
        }

        checkScans <- scanTimes[parentScan]
        rm(parentScan)

        matchPoint <- which(names(allMasses) %in% names(checkScans))
        bounds <- (matchPoint[which.min(matchPoint)] - scanInterval):
            (matchPoint[which.max(matchPoint)] + scanInterval)
        checkMasses <- allMasses[bounds]
        checkIntensities <- allIntensities[bounds]

        purityTable <- list()
        for(i in seq_along(checkMasses)) {

            curScanMzs <- checkMasses[[i]]
            curScanIntensities <- checkIntensities[[i]]
            parentIntensity <- parentIntensities[curMs2]

            mzIndex <- abs(curScanMzs - curMz)/curMz * 10^6 < ppm

            ## This will skip this step if there are no mz matches in this scan
            if(!any(mzIndex)) {
                next
            }

            checkPeaks <- curScanMzs > (curScanMzs[mzIndex] - isoWindow) &
                curScanMzs < (curScanMzs[mzIndex] + isoWindow)

            windowData <- data.frame(mzs = curScanMzs[checkPeaks],
                                   intensity = curScanIntensities[checkPeaks])

            windowData <- windowData[order(windowData$intensity, decreasing = T),] %>%
                dplyr::mutate(purity = intensity/sum(intensity)) %>%
                dplyr::filter(purity > 0.01)


            ## This function would remove C13 peak from the parent ion in the data
            ## and it would integrate the purity values.

            windowData <- check13c(windowData, ppm)
            windowData <- windowData[which.min(abs(windowData$mzs - curMz)/curMz * 10^6),]
            scan <- scanTimes[bounds[i]]

            # each table represents degree of contamination of isolation window
            # based on the other things within the scan
            purityTable[[i]] <- suppressWarnings(cbind(windowData,
                                                       parentID = parentMatches[curMs2],
                                                       scanTime = scanTimes[bounds[i]]))

        }

        purityTable <- purityTable[!sapply(purityTable, is.null)] %>%
            Reduce(rbind, .)

        scansMassSpectra[[curMs2]] <- purityTable

    }


    return(scansMassSpectra)
}
