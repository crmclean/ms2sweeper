groupSpectra <- function(groupedSpectra, massCol, mzDiff) {

    assertthat::is.string(massCol)
    assertthat::is.number(mzDiff)

    mzCol <- grep(massCol, colnames(groupedSpectra))
    if(length(mzCol) == 0) {
        stop("massCol string was not matched by groupedSpectra column names.")
    }

    groupedSpectra <- groupedSpectra[order(groupedSpectra[,mzCol]),]

    selectCooccuringPeaks <- rle(diff(groupedSpectra[,mzCol]) < mzDiff)
    groups <- vector("numeric")
    start <- 1
    groupedSpectra$groupIndex <- vector("numeric", length = nrow(groupedSpectra))
    groupCounter <- 1
    for(i in 1:length(selectCooccuringPeaks$lengths)) {

        if(isTRUE(selectCooccuringPeaks$values[i])) {
            end <- start + selectCooccuringPeaks$lengths[i]
            groupedSpectra$groupIndex[start:end] <- groupCounter
            groupCounter <- groupCounter + 1
            groups <- c(groups, start:end)
            start <- end
        } else {
            start <- start + selectCooccuringPeaks$lengths[i]
        }

    }
    return(groupedSpectra)
}


groupDaughters <- function(groupedSpectra) {

    selectCooccuringPeaks <- rle(diff(groupedSpectra$mzMean) < 0.5)
    groups <- vector("numeric")
    start <- 1
    groupedSpectra$groupIndex <- vector("numeric", length = nrow(groupedSpectra))
    groupCounter <- 1
    for(i in 1:length(selectCooccuringPeaks$lengths)) {

        if(isTRUE(selectCooccuringPeaks$values[i])) {
            end <- start + selectCooccuringPeaks$lengths[i]
            groupedSpectra$groupIndex[start:end] <- groupCounter
            groupCounter <- groupCounter + 1
            groups <- c(groups, start:end)
            start <- end
        } else {
            start <- start + selectCooccuringPeaks$lengths[i]
        }

    }
    return(groupedSpectra)
}

returnCurScan <- function(curScan, i) {

    scanData <- curScan %>% dplyr::mutate(mzMean = .data$ms2mz,
                                    mzSd = 0,
                                    intMean = .data$ms2int,
                                    intSd = 0,
                                    rtMean = .data$ms2rt,
                                    rtMin = .data$ms2rt,
                                    rtMax = .data$ms2rt,
                                    slope = 0,
                                    size = 1,
                                    firstScan = i,
                                    lastScan = i,
                                    sample = .data$ms2Samp,
                                    ms2Samp = NULL,
                                    ms2mz = NULL,
                                    ms2int = NULL,
                                    ms2rt = NULL,
                                   id = NULL)

    return(scanData)
}

returnScanDaughter <- function(curScan) {
    scanData <- curScan %>%
        dplyr::mutate(parent = .data$id,
               mzMeanDau = .data$mzMean,
               mzSdDau = 0,
               mzSdErr = 0,
               intMeanDau = .data$intMean,
               intSdDau = 0,
               intErrSd = 0,
               rtMean = .data$rtMean,
               rtMin = .data$rtMin,
               rtMax = .data$rtMax,
               slopeMean = .data$slope,
               slopeSd = 0,
               size = .data$size,
               mzMean = NULL,
               mzSd = NULL,
               intMean = NULL,
               intSd = NULL,
               id = NULL,
               firstScan = NULL,
               lastScan = NULL,
               groupIndex = NULL,
               family = NULL,
               slope = NULL,
               sample = NULL)

    return(scanData)
}

cos.sim <- function(A, B) {
    return( sum(A*B)/sqrt(sum(A^2)*sum(B^2)) )
}
