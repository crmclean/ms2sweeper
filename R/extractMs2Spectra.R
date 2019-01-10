#' @title extractMs2Spectra
#'
#' @param parentMatches - A numeric vector representing the indicies within
#' sampleData data.frame that represent a parent ion of interest. This should
#' have come from running the filterFragments function.
#' @param sampleData - A data.frame of all features belonging to a single
#' sample from which I want to extract specific information from.
#'
#' @description This function is designed to return all available MS2 peaks for
#' the given parent ions. The claim here is that each row represents a unique
#' MS2 element for a common feature.
#'
#' @return A list of MS2s that belong to the matched feature within a sample.
extractMs2Spectra <- function(parentMatches, sampleData) {

    ms2obs <- ms2Info(sampleData)
    ms2Tables <- list()

    for(i in seq_along(parentMatches)) {

        curParentPeak <- parentMatches[i]
        rowMatches <- which(ms2obs$ms2parentPeaks == curParentPeak)

        ## Do somethin here to handle things without matches
        if(length(rowMatches) == 0) {

            next

        } else {

            ms2Tables[[i]] <- data.frame(ms2mz = ms2obs$ms2mz[rowMatches],
                                            ms2int = ms2obs$ms2intensity[rowMatches],
                                            ms2rt = ms2obs$ms2rt[rowMatches],
                                            ms2Samp = ms2obs$ms2sample[rowMatches])
            names(ms2Tables)[i] <- parentMatches[i]

        }

    }

    return(ms2Tables)
}
