#' @title summarizeMs2s
#'
#' @param sweeperObj The ms2 sweeper object containing all relevant info.
#'
#' @description This function is designed to return some stats on the quality
#' of MS2. The idea here is to get a feel for what might be an appropriate
#' filter to clean up the data.
#'
#' @export
summarizeMs2s <- function(sweeperObj) {

    pureMs2s <- getMs2Pure(sweeperObj)

    if(length(pureMs2s) == 0) {
        stop("The pureMS2 slot is empty. Are you sure that MS2 spectra was harvested?")
    }

    peakScores <- sizes <- ms2ID <- scoreSd <- sizeSd <- vector()
    scanCount <- 1
    for(i in seq_along(pureMs2s)) {
        curMS2 <- pureMs2s[[i]]
        if(is.data.frame(curMS2)) {
            peakScores[scanCount] <- mean(curMS2$peakScore)
            sizes[scanCount] <- mean(curMS2$size)
            sizeSd[scanCount] <- sd(curMS2$size)
            ms2ID[scanCount] <- scanCount
            scoreSd[scanCount] <- sd(curMS2$peakScore)
            scanCount <- scanCount + 1
        } else {
            next
        }
    }
    summaryDf <- data.frame(peakScores, scoreSd, sizes, sizeSd, ms2ID)
    return(summaryDf)
}
