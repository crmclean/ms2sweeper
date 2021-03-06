% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extractMs2Spectra.R
\name{extractMs2Spectra}
\alias{extractMs2Spectra}
\title{extractMs2Spectra}
\usage{
extractMs2Spectra(parentMatches, sampleData)
}
\arguments{
\item{parentMatches}{- A numeric vector representing the indicies within
sampleData data.frame that represent a parent ion of interest. This should
have come from running the filterFragments function.}

\item{sampleData}{- A data.frame of all features belonging to a single
sample from which I want to extract specific information from.}
}
\value{
A list of MS2s that belong to the matched feature within a sample.
}
\description{
This function is designed to return all available MS2 peaks for
the given parent ions. The claim here is that each row represents a unique
MS2 element for a common feature.
}
