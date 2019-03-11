#' @title pureMs2
#'
#' @description This function is designed to return a high confidence MS2 table.
#' It mataches features that are found within two tables.
#'
#' @param ms2Matches - a list of all ms2 matches that belong to features
#' within the error of the input curMz and curRt values.
#' @param mzDiff - absolute mz error required match spectra across scans.
#'
#' @importFrom dplyr "%>%"
#'
#' @return a list of tables making up the linked scans within a sample.
pureMS2 <- function(ms2Matches, mzDiff) {

    ### some storage device
    linkedScans <- list()
    curScan <- nextScan <- data.frame()


    # Pairwise grouping of Scans ----------------------------------------------
    i <- 1
    nextExist <- T
    uniqueMS2 <- 1
    retainedScans <- data.frame()
    while(nextExist) {

        # initializing scans, and checking for next -------------------------------
        curScan <- ms2Matches[[i]]
        curScan$id <- i
        curScan$parent <- names(ms2Matches)[i]
        if(length(ms2Matches) == i) {

            curScan$groupIndex <- 0
            curScan$nextCheck <- FALSE
            linkedScans[[uniqueMS2]] <- curScan
            nextExist <- F
            next

        } else {

            nextScan <- ms2Matches[[i + 1]]
            nextScan$parent <- names(ms2Matches)[i+ 1]
            nextScan$id <- i + 1

        }


        # Checking adj scan similarity --------------------------------------------
        # grouping together the data
        groupedSpectra <- rbind(curScan, nextScan) %>%
            groupSpectra(massCol = "ms2mz", mzDiff)

        ## two conditions that would exit the function
        uniquePercent <- split(x = groupedSpectra, f = groupedSpectra$id) %>%
            sapply(function(x) {sum(x$ms2int[x$groupIndex != 0])/sum(x$ms2int)})
        if(max(uniquePercent) < .5) {

            groupedSpectra$nextCheck <- FALSE
            linkedScans[[uniqueMS2]] <- groupedSpectra
            uniqueMS2 <- uniqueMS2 + 1
            i <- i + 1
            next
        }
        rm(uniquePercent)

        combinedSpectra <- groupedSpectra[groupedSpectra$groupIndex > 0,]
        mzsInScans <- split(x = combinedSpectra$ms2mz, f = combinedSpectra$id)

        if(cos.sim(mzsInScans[[1]], mzsInScans[[2]]) < 0.9) {

            groupedSpectra$nextCheck <- FALSE
            linkedScans[[uniqueMS2]] <- groupedSpectra
            uniqueMS2 <- uniqueMS2 + 1
            i <- i + 1
            next

        } else {

            groupedSpectra$nextCheck <- TRUE
            linkedScans[[uniqueMS2]] <- groupedSpectra
            uniqueMS2 <- uniqueMS2 + 1
            i <- i + 1

        }

    }


    # Combining 2+ scans ------------------------------------------------------
    checkNext <- linkedScans %>% sapply(function(x) {all(x$nextCheck)})
    combinedScans <- list()
    counter <- 1
    if(any(checkNext)) {

        # grouping spectra with multiple scans
        skipCheck <- checkNext
        for(i in which(checkNext)) {

            ## skipping things that have already been added together
            if(!skipCheck[i]) {
                next
            }

            ## merging spectra that belong to a similar feature
            nextUp <- i + 1
            groupedScans <- rbind(linkedScans[[i]], linkedScans[[nextUp]])
            while(all(groupedScans$nextCheck)) {
                skipCheck[nextUp] <- FALSE
                nextUp <- nextUp + 1
                groupedScans <- rbind(groupedScans,
                                      linkedScans[[nextUp]])
            }

            groupedScans$nextCheck <- groupedScans$groupIndex <- NULL
            combinedScans[[counter]] <- unique(groupedScans)
            names(combinedScans)[counter] <- min(groupedScans$id)
            counter <- counter + 1
        }
        rm(skipCheck)

    }

    # interting unique spectra
    for(i in which(!checkNext)) {
        if(i == 1 || !checkNext[i - 1]) {
            singleScan <- linkedScans[[i]]
            singleScan$nextCheck <- singleScan$groupIndex <- NULL
            minScan <- min(singleScan$parent)
            singleScan <- singleScan[singleScan$parent == minScan,]
            combinedScans[[counter]] <- singleScan
            names(combinedScans)[counter] <- min(singleScan$id)
            counter <- counter + 1
        } else {
            next
        }
    }

    scanOrder <- order(names(combinedScans))
    linkedScans <- combinedScans[scanOrder]
    rm(combinedScans, scanOrder)

    # Returning Clean scans ---------------------------------------------------
    ## I need some column connecting the scans to the parent
    matchScanCount <- linkedScans %>% sapply(function(x) {length(unique(x$id))})
    for(i in seq_along(matchScanCount)) {

        if(matchScanCount[i] > 1) {

            curMergedMS2 <- linkedScans[[i]]
            groupedSpectra <- curMergedMS2 %>%
                groupSpectra(massCol = "ms2mz", mzDiff)

            multObs <- groupedSpectra[groupedSpectra$groupIndex >0,]
            noMatch <- groupedSpectra[groupedSpectra$groupIndex == 0,]
            mergedData <- multObs %>% dplyr::group_by(.data$groupIndex) %>%
                dplyr::summarize(parent = paste(unique(.data$parent), collapse = ";"),
                          mzMean = mean(.data$ms2mz),
                          mzSd = sd(.data$ms2mz),
                          intMean = mean(.data$ms2int),
                          intSd = sd(.data$ms2int),
                          rtMean = mean(.data$ms2rt),
                          rtMin = min(.data$ms2rt),
                          rtMax = max(.data$ms2rt),
                          slope = stats::coef(stats::lm(.data$ms2int~.data$ms2rt))[2],
                          size = length(.data$ms2int),
                          firstScan = min(.data$id),
                          lastScan = max(.data$id),
                          sample = unique(.data$ms2Samp)) %>%
                dplyr::mutate(groupIndex = NULL) %>%
                data.frame()

            noMatch <- returnCurScan(noMatch,noMatch$id) %>%
                dplyr::mutate(groupIndex = NULL)
            linkedScans[[i]] <- rbind(mergedData, noMatch)

        } else {

            linkedScans[[i]] <- returnCurScan(linkedScans[[i]],i)
        }

    }

    return(linkedScans)
}




