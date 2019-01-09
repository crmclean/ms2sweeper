#' @title sweeperObj
#'
#' @description This object is a generic object designed to run the different
#' functions of the ms2sweeper package. The slots represent content or data
#' that the package uses throughout the different functions.
#'
#' @slot features - A table containing an m/z and rt column used to select
#' specific features from ms2 data.
#' @slot ms2Path - A string vector of paths to the ms2 raw data files.
#' @slot ms2Data - The read in ms2 data.
#' @slot ms2Harvest - Harvested ms2 data for all inputed parent ions from
#' features.
#' @slot ms2Pure - list of purified ms2s based on the cleaning criterion.
setClass(
    Class = "sweeperObj",

    representation = list(

        features = "data.frame",
        ms2Path = "character",
        ms2Data = "list",
        harvestDaughters = "list",
        harvestParents = "list",
        outputPath = "character",
        ms2Pure = "list"

    )
)

setMethod(

    f="initialize",
    signature="sweeperObj",
    definition=function(.Object, features, ms2Path) {

        cat("~~~ sweeperObj: Initializator ~~~ \n")
        # Adding featureAttr obj -----------------------------------------------
        if(missing("features") | missing("ms2Path")) {
            stop("Check to make sure features and ms2Path arguement is given.")
        }


        .Object <- setFreatures(sweeperObj = .Object, features = features)
        .Object <- setMs2Path(sweeperObj = .Object, ms2Path = ms2Path)

        cat("~~~ Loading fragment data into R ~~~ \n")
        .Object <- setMs2Data(sweeperObj = .Object)

        cat("~~~ The fragments are loaded ~~~ \n")
        return(.Object) # return of the object
    }
)

createSweeper <- function(features, ms2Path) {
    sweeperObj <- new(Class="sweeperObj", features, ms2Path)
    return(sweeperObj)
}

#' @title get_set
#'
#' @description This file contains the getter and setter functions used to
#' access slots from sweeperObj object.
#'
#' @export
# Get Functions -----------------------------------------------------------
getFreatures <- function(sweeperObj) {
    return(sweeperObj@features)
}

getMs2Path <- function(sweeperObj) {
    return(sweeperObj@ms2Path)
}

getMs2Data <- function(sweeperObj) {
    return(sweeperObj@ms2Data)
}

getHarvestDaughters <- function(sweeperObj) {
    return(sweeperObj@harvestDaughters)
}

getHarvestParents <- function(sweeperObj) {
    sweeperObj@harvestParents
}

getMs2OutputPath <- function(sweeperObj) {
    sweeperObj@outputPath
}

getMs2Pure <- function(sweeperObj) {
    sweeperObj@ms2Pure
}


# Set Functions -----------------------------------------------------------
setFreatures <- function(sweeperObj, features) {

    if(!is.data.frame(features)) {
        stop("Make sure input to features arguement is a data.frame")
    }
    if(!all(c("mz", "rt") %in% colnames(features))) {
        stop("Make sure input to features arguement has columns named \n
             mz and rt.")
    }
    if(!(is.numeric(features$mz) & is.numeric(features$rt))) {
        stop("Make sure mz and rt columns in features arguement are numeric.")
    }

    sweeperObj@features <- features
    return(sweeperObj)
    }

setMs2Path <- function(sweeperObj, ms2Path) {

    if(!all(file.exists(ms2Path))) {
        stop("Make sure ms2Path arguement points to real file.")
    }

    sweeperObj@ms2Path <- ms2Path
    return(sweeperObj)
}

setMs2Data <- function(sweeperObj) {

    fragmentation_list <- list()
    ms2Files <- getMs2Path(sweeperObj)
    length(fragmentation_list) <- length(ms2Files)
    for(i in 1:length(ms2Files)) {

        suppressMessages(theFragments <- suppressMessages(xcms::xcmsSet(ms2Files[i],
                                                                        method="MS1") %>%
                                                              xcms::xcmsFragments(., snthresh = 1)))

        theFragments <- data.frame(theFragments@peaks) %>%
            dplyr::mutate(dataFile = ms2Files[i])

        fragmentation_list[[i]] <- theFragments
    }
    sweeperObj@ms2Data <- fragmentation_list
    return(sweeperObj)
}

setHarvestDaughters <- function(sweeperObj, harvestDaughters) {
    if(is.list(harvestDaughters)) {
        sweeperObj@harvestDaughters <- harvestDaughters
    } else {
        stop("Incorrect input for harvestDaughter Arg w/in HarvestMs2 Function.")
    }
    return(sweeperObj)
}

setHarvestParents <- function(sweeperObj, harvestParent) {
    if(is.list(harvestParent)) {
        sweeperObj@harvestParents <- harvestParent
    } else {
        stop("Incorrect input for harvestParent Arg w/in HarvestMs2 Function.")
    }
    return(sweeperObj)
}

setMs2Pure <- function(sweeperObj, ms2List) {

    if(is.list(ms2List)) {
        sweeperObj@ms2Pure <- ms2List
    } else {
        stop("Incorrect input for ms2Pure slot.")
    }
    return(sweeperObj)
}

# Consider adding a clear MS2 function ------------------------------------


clearMs2Data <- function(sweeperObj) {
    sweeperObj@ms2Data <- list()
    return(sweeperObj)
}
