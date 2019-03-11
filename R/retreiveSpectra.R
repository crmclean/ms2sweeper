#' @title retreiveSpectra
#'
#' @description This function is desgined to output MS2 spectra gathered from
#' the sweeper algorithm in feature (mz/rt) specfic tables.
#'
#' @param sweeperObj - sweeper object that has been used throughout the
#' algorithm.
#' @param outputPath - path to store output files to run on metfrag.
#'
#' @export
retreiveSpectra <- function(sweeperObj, outputPath) {

    assertthat::is.dir(outputPath)

    curWd <- getwd()
    setwd(outputPath)

    ms2Scans <- getMs2Pure(sweeperObj)

    if(length(ms2Scans) == 0) {
        stop("There is no pure MS2 scans within sweeperObj to output.")
    }

    for(i in seq_along(ms2Scans)) {

        curMS2 <- ms2Scans[[i]]

        if(!is.data.frame(curMS2)) {
            next
        }

        if(nrow(curMS2) > 0) {

            parMz <- signif(curMS2$parMz[1])
            parRt <- signif(curMS2$parRt[1])
            outTable <- curMS2[,grep("mzMeanDau|intMeanDau", colnames(curMS2))]
            colnames(outTable) <- c("mz", "intensity")
            utils::write.table(paste0(parMz, "_", parRt, "_ms2.txt"),
                        x = outTable,
                        row.names = F,
                        col.names = F)

        }


    }

    return(NULL)

}
