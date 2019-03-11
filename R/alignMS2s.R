#' @title alignMS2s
#'
#' @description This function is designed to align MS2 spectra from a common
#' parent ion. The input consists of the sweeperObj and out output is the
#' same object with the pureMS2 slot filled.
#'
#' @param sweeperObj - sweeper object containing ms2 data for each matched
#' feature per sample
#' @param ppm - ppm error to check for similarities between parent ions across
#' samples.
#' @param betweenSampThresh - minimum p probability value required to retain a
#' peak during across sample comparison.
#' @param poisThresh - minimum p value required to retaion a peak during
#' modeling of occurance frequencies.
#' @param cores - number of cores running the algorithm.
#' @param mzDiff - absolute mass difference required to group common peaks of
#' across distinct ms2 scans.
#'
#' @importFrom foreach "%dopar%"
#' @importFrom dplyr "%>%"
#' @importFrom stats "sd"
#'
#' @return the sweeper object with the pureMS2 slot filled.
#' @export
alignMS2s <- function(sweeperObj, ppm = 5, betweenSampThresh = 0.15,
                      poisThresh = 0.05, cores = 1, mzDiff = 0.5) {


    geomean <- function(x, na.rm=TRUE) {
        exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
    }

    parents <- getHarvestParents(sweeperObj)
    daughters <- getHarvestDaughters(sweeperObj)
    names(daughters) <- daughters %>% lapply(function(x) {as.character(x$sample[1])})

    if(length(parents) != length(daughters)) {
        stop("Parents do not equal daughters. Check harvestMS2s function.")
    }

    allParents <- Reduce(rbind, parents) %>% .[order(.$mzs),]
    selectCooccuringPeaks <- rle(diff(allParents$mzs)/allParents$mzs[2:nrow(allParents)]*10^6 < ppm)
    groups <- vector("numeric")
    start <- 1
    allParents$groupIndex <- vector("numeric", length = nrow(allParents))
    groupCounter <- 1
    for(i in 1:length(selectCooccuringPeaks$lengths)) {

        if(isTRUE(selectCooccuringPeaks$values[i])) {
            end <- start + selectCooccuringPeaks$lengths[i]
            allParents$groupIndex[start:end] <- groupCounter
            groupCounter <- groupCounter + 1
            groups <- c(groups, start:end)
            start <- end
        } else {
            start <- start + selectCooccuringPeaks$lengths[i]
        }

    }

    ## this part of the algo could be written in parallel
    xx <- allParents %>% split(f = allParents$groupIndex)


    # Fixing things that only occour once -------------------------------------
    orig <- groupCounter
    zeroList <- xx[[1]] %>% split(f = paste(xx[[1]]$family, xx[[1]]$sample))
    for(i in seq_along(zeroList)) {

        zeroList[[i]]$groupIndex <- groupCounter
        groupCounter <- groupCounter + 1
    }

    xx[orig:(groupCounter-1)] <- zeroList
    xx <- xx[-1]
    rm(zeroList, orig)

    doMC::registerDoMC(cores)
    allMS2s <- foreach::foreach(i = seq_along(xx)) %dopar% {

        message(i)
        storeMS2s <- list()
        ms2Counter <- 1

        ## solving this for i first, then doing the rest later
        curGroup <- xx[[i]]
        curGroup$scanId <- paste(curGroup$family,
                                 as.character(curGroup$sample),
                                 sep = "_")

        sampleGroups <- curGroup %>% split(f = curGroup$scanId)
        minScanTime <- sampleGroups %>% sapply(function(x) {min(x$rt)})
        sampleGroups <- sampleGroups[order(minScanTime)]
        rm(minScanTime)

        # extracting relevant ms2 data
        ms2Matches <- list()
        for(j in seq_along(sampleGroups)) {
            checkDaughter <- daughters[[which(names(daughters) %in% sampleGroups[[j]]$sample[1])]]
            ms2Matches[[j]] <- checkDaughter[checkDaughter$family == sampleGroups[[j]]$family[1],]
        }
        rm(j)

        ## merging ms2 data across samples
        concensusMS2s <- pureDaughters(ms2Matches, mzDiff)

        # Bayesian prob test ------------------------------------------------------
        for(j in seq_along(concensusMS2s)) {

            curMS2 <- concensusMS2s[[j]]

            if(any(grepl(curMS2$parent, pattern = ";"))) {
                curPars <- unique(unlist(strsplit(curMS2$parent, ";")))
                curPars <- as.numeric(curPars)
            } else {
                curPars <- as.numeric(unique(curMS2$parent))
            }


            parPurity <- sampleGroups[curPars] %>%
                sapply(function(x) {mean(x$purity)})

            parIntensity <- sampleGroups[curPars] %>%
                sapply(function(x) {mean(x$intensity)})

            parMz <- sampleGroups[curPars] %>%
                sapply(function(x) {mean(x$mzs)})

            parRt <- sampleGroups[curPars] %>%
                sapply(function(x) {mean(x$rt)})

            parProbs <- parIntensity*parPurity/sum(parIntensity * parPurity)

            curMS2$probabilty <- 0
            curMS2$peakScore <- 0
            geoMeanInts <- curMS2$intMeanDau/geomean(curMS2$intMeanDau)
            peakProbs <- 1/(1+1/geoMeanInts)

            ## doing bayesian stuff here
            for(kk in 1:nrow(curMS2)) {

                curRow <- curMS2$parent[kk]
                curPeakProb <- peakProbs[kk]
                if(any(grepl(";", curRow))) {
                    peakOrigs <- as.numeric(unlist(strsplit(curRow, ";")))
                } else {
                    peakOrigs <- as.numeric(curRow)
                }

                checkPars <- which(curPars %in% peakOrigs)

                curProb <- sum(parProbs[checkPars])
                curMS2$probabilty[kk] <- curProb*curPeakProb/(
                    curProb*curPeakProb + (1-curProb)*(1-curPeakProb))
                curMS2$peakScore[kk] <- sum((parIntensity*parPurity)[checkPars]) *
                    curMS2$size[kk]

            }

            cleanMs2 <- curMS2[curMS2$probabilty > betweenSampThresh,]

            ## doing poisson stuff here
            if(sum(cleanMs2$size > 1) > 1) {
                expectedValue <- mean(cleanMs2$size[cleanMs2$size > 1])
                keepPois <- ppois(q = cleanMs2$size, lambda = expectedValue) > poisThresh
                cleanMs2 <- cleanMs2[keepPois,]
            }

            cleanMs2$parMz <- mean(parMz)
            cleanMs2$parMzSd <- sd(parMz)
            cleanMs2$parRt <- mean(parRt)
            cleanMs2$parRtSd <- sd(parRt)

            if(nrow(cleanMs2) > 0) {
                storeMS2s[[ms2Counter]] <- cleanMs2
                ms2Counter <- ms2Counter + 1
            }

        }

        return(storeMS2s)

    }

    allMS2s <- allMS2s %>% unlist(recursive=F)
    sweeperObj <- setMs2Pure(sweeperObj = sweeperObj, ms2List = allMS2s)

    return(sweeperObj)

}
