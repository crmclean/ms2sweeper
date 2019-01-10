#' @title pureDaughters
#'
#' @description This function works the same as the pureMS2 function, but is
#' slightly different for matching spectra across samples. This only means
#' that it has a few more error messages.
pureDaughters <- function(ms2Matches, mzDiff) {

    # Merging Scans by Similarity ---------------------------------------------
    linkedScans <- list()
    curScan <- data.frame()
    i <- 1
    nextExist <- T
    uniqueMS2 <- 1
    scanGroups <- list()
    while(nextExist) {

        # initializing scans, and checking for next -------------------------------
        curScan <- ms2Matches[[i]]
        curScan$curGroup <- 0
        curScan$id <- i
        curScan$parent <- names(ms2Matches)[i]
        if(length(ms2Matches) == i) {

            curScan$curGroup <- uniqueMS2
            linkedScans[[uniqueMS2]] <- curScan
            nextExist <- F
            next

        } else if(i == 1) {

            curScan$curGroup <- uniqueMS2
            linkedScans[[uniqueMS2]] <- curScan
            uniqueMS2 <- uniqueMS2 + 1
            i <- i + 1
            next

        }

        ## checking scans accross groups
        nextGroup <- TRUE
        addedSpectra <- FALSE
        j <- 1
        while(nextGroup) {

            ### add another bool here to check for things that are already linked within data
            curMs2Group <- linkedScans[[j]]
            curScan$groupIndex <- NULL

            groupedSpectra <- rbind(curScan, curMs2Group) %>%
                .[order(.$mzMean),] %>%
                groupSpectra(massCol = "mzMean", mzDiff)

            curScan <- groupedSpectra[groupedSpectra$curGroup == 0,]
            curMs2Group <- groupedSpectra[groupedSpectra$curGroup > 0,]

            # finishing up first similarity check
            uniquePercent <- sum(curScan$intMean[curScan$groupIndex > 0]/sum(curScan$intMean))
            if(uniquePercent < 0.5) {

                if(j == length(linkedScans)) {

                    if(!addedSpectra) {
                        curScan$curGroup <- uniqueMS2
                        curScan$groupIndex <- NULL
                        linkedScans[[uniqueMS2]] <- curScan
                        uniqueMS2 <- uniqueMS2 + 1
                    }
                    i <- i + 1
                    nextGroup <- FALSE

                } else {
                    j <- j + 1
                }
                next
            }

            checkGroups <- curScan$groupIndex[curScan$groupIndex != 0]
            newGroupedMzs <- curScan$mzMean[curScan$groupIndex %in% checkGroups]


            combinedSpectra <- curMs2Group[curMs2Group$groupIndex %in% checkGroups,]
            groupedMzs <- combinedSpectra %>%
                group_by(groupIndex) %>%
                summarise(mzMean = mean(mzMean)) %>%
                data.frame()
            groupedMzs <- groupedMzs$mzMean

            if(length(groupedMzs) < length(newGroupedMzs)) {

                checkCurScan <- checkGroups[duplicated(checkGroups)]

                # fixing samples where multiple peaks belong to the same group
                for(k in seq_along(checkCurScan)) {

                    curMasses <- curScan$mzMean[curScan$groupIndex %in% checkCurScan[k]]
                    oldMasses <- curMs2Group$mzMean[curMs2Group$groupIndex %in% checkCurScan[k]]
                    minDistIndex <- which.min(abs(curMasses - mean(oldMasses)))
                    curScan$groupIndex[curScan$groupIndex %in% checkCurScan[k]][-minDistIndex] <- 0

                }

                #combinedSpectra$[combinedSpectra$groupIndex %in% checkCurScan,]
                warning(paste("mzDiff may be too large.",
                              "You are grouping peaks that belong to the",
                              "same sample.\n Consider lowering mzDiff param."))

                newGroupedMzs <- curScan$mzMean[curScan$groupIndex %in% checkGroups]
            }

            cosScore <- cos.sim(newGroupedMzs, groupedMzs)

            # checking similarity in peak intensity
            if(cosScore < 0.9) {

                if(j == length(linkedScans)) {

                    if(!addedSpectra) {
                        curScan$curGroup <- uniqueMS2
                        linkedScans[[uniqueMS2]] <- curScan
                        uniqueMS2 <- uniqueMS2 + 1
                    }
                    i <- i + 1
                    nextGroup <- FALSE

                } else {
                    j <- j + 1
                }
                next
            }

            # adding old spectra to old features
            addedSpectra <- TRUE
            groupedSpectra$groupIndex <- NULL
            groupedSpectra$curGroup <- max(groupedSpectra$curGroup)
            linkedScans[[j]] <- groupedSpectra

            if(j < length(linkedScans)) {
                j <- j + 1
            } else {
                i <- i + 1
                nextGroup <- FALSE
            }

        }
    }

    # Returning Clean scans ---------------------------------------------------
    ## I need some column connecting the scans to the parent
    matchScanCount <- linkedScans %>% sapply(function(x) {length(unique(x$id))})
    for(i in seq_along(matchScanCount)) {


        ## debug some error here...
        if(matchScanCount[i] > 1) {

            groupedSpectra <- groupSpectra(linkedScans[[i]],
                                           massCol = "mzMean",
                                           mzDiff)

            groupedSpectra$curGroup <- NULL
            multObs <- groupedSpectra[groupedSpectra$groupIndex >0, ]
            noMatch <- groupedSpectra[groupedSpectra$groupIndex == 0,]
            mergedData <- multObs %>%
                dplyr::group_by(groupIndex) %>%
                dplyr::summarize(parent = paste(unique(id), collapse = ";"),
                          mzMeanDau = mean(mzMean),
                          mzSdDau = sd(mzMean),
                          mzSdErr = sd(mzSd),
                          intMeanDau = mean(intMean),
                          intSdDau = sd(intMean),
                          intErrSd = sd(intSd),
                          rtMean = mean(rtMean),
                          rtMin = min(rtMin),
                          rtMax = max(rtMax),
                          slopeMean = mean(slope),
                          slopeSd = sd(slope),
                          size = sum(size)) %>%
                dplyr::mutate(groupIndex = NULL) %>%
                data.frame()

            noMatch <- returnScanDaughter(noMatch)
            linkedScans[[i]] <- rbind(mergedData, noMatch)

        } else {

            linkedScans[[i]] <- returnScanDaughter(linkedScans[[i]])
        }

    }

    return(linkedScans)
}




