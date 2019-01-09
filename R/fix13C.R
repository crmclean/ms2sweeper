#' @title fix13C
#'
#' @description This function is designed to remove peaks that may be from 13C
#' peak within the MS2 scan
#'
#' @param ms2Table - ms2 spectra table being checked for 13C peaks.
#'
#' @return A ms2 spectra table corrected for 13C peaks.
#' @export
fix13C <- function(ms2Table) {

    ms2Table$index <- 1:nrow(ms2Table)

    ## checking for presence of c13 peaks
    peakDistance <- dist(ms2Table$ms2mz) %>%
        as.matrix()
    peakDistance[lower.tri(peakDistance, diag = TRUE)] <- 0
    peakDistance <- subset(reshape2::melt(peakDistance), value > 0) %>%
        .[order(.$Var1),]

    ### I SHOULD CHANGE THIS LATER TO BE AN IMPUT VALUE
    pairCheck <- peakDistance$value < 1.3

    ## removing c13 peak
    if(any(pairCheck)) {

        c13match <- which(pairCheck)
        # goes through and corrects each pair
        for(pair in seq_along(c13match)) {
            feature1 <- peakDistance$Var1[c13match[pair]]
            feature2 <- peakDistance$Var2[c13match[pair]]

            # checking to make sure a feature hasn't already been removed
            if(!all(c(feature1, feature2) %in% ms2Table$index)) {
                next
            }

            c13 <- which.min(c(ms2Table$ms2int[feature1], ms2Table$ms2int[feature2]))

            if(c13 == 1) {
                ms2Table <- ms2Table[-feature1,]
            } else {
                ms2Table <- ms2Table[-feature2,]
            }
        }

    }

    ms2Table$index <- NULL
    return(ms2Table)

}

