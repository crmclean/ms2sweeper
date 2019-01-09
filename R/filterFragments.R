#' @title filterFragments
#'
#' @param curMz - current m/z value of feature searched within the data
#' @param curRt - current rt value of feature searched within the data
#' @param sampleData - current sample being screened for ms2 matches
#' @param ppm - mass error match between sample m/zs and entered curMz
#' @param rtError - rt error between sample feature's rt and curRt
#'
#' @description This function is designed to find mz and rt matches between the
#' current feature being investigated and a sample. It shoud return indicies of
#' features that are within the allowed error.
#'
#' @return a vector of matching indicies.
#' @export
filterFragments <- function(curMz, curRt, sampleData, ppm, rtError) {


    # Matching things by ppm error --------------------------------------------
    # calculating ppm
    matchingRowsMS2s <- which(abs(sampleData$mz - curMz)/curMz * 10^6 < ppm)
    parentPeaks <- sampleData$msLevel[matchingRowsMS2s] == 1

    ## figure out what to return so it interacts with the past section
    if(!any(parentPeaks)) {
        return(NULL)
    }

    matchingRowsMS2s <- matchingRowsMS2s[parentPeaks]


    # Matching things by rt error ---------------------------------------------
    ## checking rt error window
    matchedParents <- matchingRowsMS2s[
        abs(curRt - sampleData$rt[matchingRowsMS2s]) < rtError
        ] %>% unique()

    if(length(matchedParents) == 0) {
        return(NULL)
    }

    ## returning a numeric vector index
    return(matchedParents)

}
