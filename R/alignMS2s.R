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
#'
#' @importFrom foreach "%dopar%"
#'
#' @return the sweeper object with the pureMS2 slot filled.
#' @export
alignMS2s <- function(sweeperObj, ppm = 5, betweenSampThesh = 0.05,
                      poisThresh = 0.05, cores = 1, mzDiff = 0.15) {

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

    doMC::registerDoMC(cores)
    allMS2s <- foreach::foreach(i = seq_along(xx)) %dopar% {

        storeMS2s <- list()
        ms2Counter <- 1

        ## solving this for i first, then doing the rest later
        curGroup <- xx[[i]]
        curGroup$scanId <- paste(curGroup$family,
                                 as.character(curGroup$sample),
                                 sep = "_")

        sampleGroups <- curGroup %>% split(f = curGroup$scanId)
        minScanTime <- sampleGroups %>% sapply(function(x) {min(x$scanTime)})
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

        ## computing properties about the parent
        parInts <- sampleGroups %>%
            sapply(function(x) {sum(x$intensity)})

        parMz <- sampleGroups %>%
            sapply(function(x) {mean(x$mzs)})

        parRt <- sampleGroups %>%
            sapply(function(x) {mean(x$scanTime)})

        meanParPur <- sampleGroups %>%
            sapply(function(x) {mean(x$purity)})

        parProbs <- parInts*meanParPur/sum(parInts * meanParPur)

        ## cleaning up features with some stats
        for(j in seq_along(concensusMS2s)) {

            curMS2 <- concensusMS2s[[j]]
            if(is.character(curMS2$parent)) {
                parents <- curMS2$parent %>% sapply(strsplit, split = ";")
                uniqueParents <- unlist(parents) %>% unique() %>% as.numeric()
                keepPeaks <- parents %>%
                    sapply(function(x) {
                        sum(parProbs[as.numeric(x)])
                    }) > betweenSampThesh

                curMS2 <- curMS2[keepPeaks,]
            } else {
                uniqueParents <- curMS2$parent %>% unique()
            }

            keepPois <- ppois(curMS2$size, mean(curMS2$size)) > poisThresh
            curMS2 <- curMS2[keepPois,]

            curMS2$parentMz <- mean(parMz[uniqueParents])
            curMS2$parentMzSd <- sd(parMz[uniqueParents])
            curMS2$parentRt <- mean(parRt[uniqueParents])
            curMS2$parentRtSd <- sd(parRt[uniqueParents])

            storeMS2s[[ms2Counter]] <- curMS2
            ms2Counter <- ms2Counter + 1
        }

        return(storeMS2s)

    }

    allMS2s <- allMS2s %>% unlist(recursive=F)
    sweeperObj <- setMs2Pure(sweeperObj = sweeperObj, ms2List = allMS2s)

    return(sweeperObj)

}
