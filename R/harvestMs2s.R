#' @title harvestMS2
#'
#' @description This function connects all MS2 scans belonging to a single feature
#' within a single data sample. It returns the sweeperObj with a slot for the
#' parent ions and a slot for the merged MS2s filled up. It also does a 13C
#' correction of the MS2s if a 13C peak is detected within the parent ion's
#' scan.
#'
#' @param sweeperObj - object containing all MS2 data that will be mined.
#' @param rtError - Time window to check for matches in rt betwee parent ions
#' from XCMS and those from the MS2 fragment table.
#' @param ppm - Error window for features.
#' @param isoWindow - isolation window for MS2s
#' @param clearData - Boolean. Deletes raw data from sweepersObj. Good to
#' reduce amount of memory taken up by sweeperObj.
#' @param mzDiff - mz error window used to group peaks from spectra together.
#' that belong to the same features but come in different scans.
#' @param cores - number of cores to use during harvest function.
#'
#' @importFrom foreach "%dopar%"
#' @importFrom dplyr "%>%"
#'
#' @return sweeperObj with slots for parents and daughters filled.
#' @export
harvestMS2 <- function(sweeperObj, ppm = 5, rtError = 15,
                       isoWindow = 2,
                       cores = 3,
                       clearData = TRUE,
                       mzDiff = 0.5) {

    if(length(getMs2Data(sweeperObj)) == 0) {
        stop("MS2 data has not been loaded into sweeper object.")
    }

    allMs2s <- getMs2Data(sweeperObj)
    features <- getFreatures(sweeperObj)

    # Looping through samples -------------------------------------------------

    cl <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)
    sampData <- foreach::foreach(i = seq_along(allMs2s)) %dopar% {

        message(paste("Currently on sample:", i))

        curSamp <- allMs2s[[i]]$dataFile[1]

        # Loading in file info ----------------------------------------------------
        rawData <- MSnbase::readMSData(files = curSamp,
                                       mode = "onDisk",
                                       msLevel. = 1)
        scanTimes <- xcms::rtime(rawData)
        allMasses <- xcms::mz(rawData)
        allIntensities <- xcms::intensity(rawData)
        rm(rawData,curSamp)

        # objects being initialized should be stored together in a list
        lastMz <- lastRt <- 0
        storeParent <- storeDaughter <- list()
        matchCounter <- family <- 1


        # Checking each mz value individually -------------------------------------
        for(j in 1:nrow(features)) {

            # checking if I need to skip ----------------------------------------------
            curMz <- features$mz[j]
            curRt <- features$rt[j]

            ## checking to see if current mz is equivalent to last
            if(identical(c(curMz, curRt), c(lastMz, lastRt))) {
                next
            }

            # filtering data by mz and rt error ---------------------------------------
            # this only checks against ms1 spectra
            parentMatches <- filterFragments(curMz, curRt,
                                             sampleData = allMs2s[[i]],
                                             ppm, rtError)

            # this will be null if no match between ms2 data and parent ion is
            # made
            if(is.null(parentMatches)) {
                next
            }

            # extracting MS2 spectra for each peak --------------------------------------
            ms2Matches <- extractMs2Spectra(parentMatches = parentMatches,
                                            sampleData = allMs2s[[i]])

            parentMatches <- parentMatches[which(parentMatches %in% names(ms2Matches))]


            # this will be null if no ms2 peaks were found for a given peak
            nullMatch <- sapply(ms2Matches, is.null)
            if(all(nullMatch)) {
                next
            }

            ms2Matches <- ms2Matches[!nullMatch]

            # checking rt purity ------------------------------------------------------
            purityTable <- parentPurity(sampleData = allMs2s[[i]],
                                        curMz = curMz,
                                        parentMatches = parentMatches,
                                        allMasses,
                                        allIntensities,
                                        scanTimes,
                                        ppm,
                                        isoWindow)

            # this will be null if...
            if(is.null(purityTable)) {
                next
            }

            ## cleans up MS2s for isotope peaks
            for(k in seq_along(purityTable)) {

                ms2Match <- which(names(ms2Matches) == unique(purityTable[[k]]$parentID))
                if(any(purityTable[[k]]$has13C)) {

                    if(length(ms2Match) != 1) {
                        next()
                    } else {
                        ms2Matches[[ms2Match]] <- fix13C(ms2Table = ms2Matches[[ms2Match]])
                    }

                }

                lineTermsPurity <- coef(lm(data = purityTable[[k]],
                                     formula = purity ~ scanTime))
                scanPurity <- lineTermsPurity[1] +
                    lineTermsPurity[2] * ms2Matches[[ms2Match]]$ms2rt[1]
                lineTermsIntensity <- coef(lm(data = purityTable[[k]],
                                         formula = intensity ~ scanTime))
                scanIntensity <- lineTermsIntensity[1] +
                    lineTermsIntensity[2] * ms2Matches[[ms2Match]]$ms2rt[1]
                purityTable[[k]] <- data.frame(mzs = purityTable[[k]]$mzs[1],
                           intensity = scanIntensity,
                           purity = scanPurity,
                           parentID = purityTable[[k]]$parentID[1],
                           rt = ms2Matches[[ms2Match]]$ms2rt[1])

            }

            mergedMS2s <- pureMS2(ms2Matches, mzDiff)

            ## there is an issue here... I need to work on this. Maybe get rid of family

            ## merging purity table w/ ms2 spcctra
            parents <- daughters <- list()
            for(k in seq_along(mergedMS2s)) {

                parentsIds <- mergedMS2s[[k]]$parent %>% strsplit(";") %>%
                    unlist() %>%
                    unique()
                parents[[k]] <- purityTable[parentMatches %in% parentsIds] %>%
                    Reduce(rbind,.) %>%
                    dplyr::mutate(family, has13C = NULL)
                daughters[[k]] <- mergedMS2s[[k]] %>% dplyr::mutate(family)
                family <- family + 1

            }


            storeParent[[matchCounter]] <- Reduce(rbind, parents)
            storeDaughter[[matchCounter]] <- Reduce(rbind, daughters)
            matchCounter <- matchCounter + 1

            # update general output storage -------------------------------------------
            lastMz <- curMz
            lastRt <- curRt

        }

        # Storing MS2 List --------------------------------------------------------
        tempPar <- Reduce(rbind, storeParent)
        tempDaughter <- Reduce(rbind, storeDaughter)
        tempPar$sample <- tempDaughter$sample[1]

        rm(storeDaughter, storeParent)

        return(list(tempPar, tempDaughter))
    }

    parseParents <- parseDaughters <- list()
    for(i in seq_along(sampData)) {
        parseParents[[i]] <- sampData[[i]][[1]]
        parseDaughters[[i]] <- sampData[[i]][[2]]
    }
    rm(sampData)

    sweeperObj <- setHarvestParents(sweeperObj, parseParents)
    sweeperObj <- setHarvestDaughters(sweeperObj, parseDaughters)

    if(clearData) {
        sweeperObj <- clearMs2Data(sweeperObj)
    }

    return(sweeperObj)
}
