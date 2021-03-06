% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SumColony.R
\name{SumColony}
\alias{SumColony}
\title{Summarize the results from Colony2 analysis for parentage on muliple PseudoBabies simulations}
\usage{
SumColony(
  nSim,
  Miss_data = T,
  MissingVec = c(0, 2.5),
  Geno_Error = T,
  ErrorVals = 1,
  Cutoff = 0.8,
  Markerfile = "MarkerPanels.csv"
)
}
\arguments{
\item{Miss_data}{Would you like to add missing data to the simulations? Defaults to FALSE}

\item{MissingVec}{A vector of values for the percent of missing data in the input file. Required if Miss_data=TRUE. Defaults to 0 and 2.5 percent.}

\item{Geno_Error}{Would you like to add genotyping error to the simualated data? Defaults to FALSE}

\item{ErrorVals}{How many error values were specified in Loci_error.csv? Required if Geno_Error=TRUE. Defaults to 1 percent.}

\item{Cutoff}{cutoff value for probability value accepting parentage.}

\item{Markerfile}{.csv for marker panels that will be evaluated.}

\item{nSims}{Number of replicate simulations}
}
\description{
This function summarizes the results from multiple runs of Colony2 for parentage inference from PseudoBabies simulations.
}
\examples{
SumColony()
}
\keyword{Colony}
\keyword{parentage}
