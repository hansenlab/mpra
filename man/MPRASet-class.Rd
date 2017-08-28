\name{MPRASet-class}
\Rdversion{1.1}
\docType{class}
\alias{MPRASet-class}
\alias{show,MPRASet-method}
\alias{MPRASet}
\alias{getRNA}
\alias{getDNA}
\title{Class \code{"MPRASet"}}
\description{
  A container for data from massively parallel reporter assays
  (MPRA). Builds on the SummarizedExperiment class.
}
\usage{
## Constructor
MPRASet(DNA = new("matrix"), RNA = new("matrix"),
        barcode = new("DNAStringSet"), eid = new("character"),
        eseq = new("DNAStringSet"), ...)

## Accessors
getRNA(object, aggregate = FALSE)
getDNA(object, aggregate = FALSE)
}
\arguments{
  \item{object}{A \code{MPRASet} object.}
  \item{aggregate}{A \code{logical} indicating if data should be
    aggregated to the element level (by summing across barcodes).}
  \item{DNA}{}
  \item{RNA}{}
  \item{barcode}{}
  \item{eid}{}
  \item{eseq}{}
  \item{...}{}
}
\section{Slots}{
  Slots are as described in a \code{SummarizedExperiment}.  Of
  particular interest are \code{colData} which describes the phenotype
  data.  The \code{assay} slot holds the assay data, with specific assay
  names \code{RNA} and \code{DNA} (accessed by \code{getRNA} and
  \code{getDNA}). Element and barcode data are accessible in the
  \code{rowData} slot.
}
\section{Extends}{
Class \code{"\linkS4class{SummarizedExperiment}"}, directly.
}
\section{Accessors}{
\describe{
\item{\code{getDNA}:}{Gets the DNA channel data.}
\item{\code{getRNA}:}{Gets the RNA channel data.}
}
}
\references{
  FIXME
}
\seealso{
  \code{\linkS4class{SummarizedExperiment}} for the basic class that is
  used as a building block.
}
\examples{
showClass("MPRASet")
}
\keyword{classes}
