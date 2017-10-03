\name{fit_corr}
\alias{fit_corr}
\title{Paired linear model fit for MPRA differential analysis}
\description{
Fits a mixed effects model (random intercept) for each element as part
of the differential analysis routine. This accounts for pairing between
the columns of the \code{MPRASet} object.
}
\usage{
fit_corr(object, aggregate = c("mean", "sum", "none"), design,
         normalize = TRUE, block = NULL, plot = TRUE, span = 0.4, ...)
}
\arguments{
  \item{object}{An object of class \code{MPRASet}.}
  \item{aggregate}{Aggregation method over barcodes: \code{"mean"} to
    use the average of barcode-specific log ratios, \code{"sum"} to use
    the log ratio of summed RNA and DNA counts, \code{"none"} to perform
    no aggregation (counts have already been summarized over barcodes).} 
  \item{design}{Design matrix specifying comparisons of interest.}
  \item{normalize}{If \code{TRUE}, perform total count normalization
    before model fitting.}
  \item{block}{A vector giving the sample designations of the columns of
    \code{object}. The default, \code{NULL}, indicates that all columns
    are separate samples.} 
  \item{plot}{If \code{TRUE}, plot the copy number-variance relationship.}
  \item{span}{The smoothing span for \code{lowess} in estimating the
    copy number-variance relationship. Default: 0.4.} 
  \item{\dots}{Further arguments to be passed to \code{lmFit} for
    obtaining residual standard deviations used in estimating the
    copy number-variance relationship.}
}
\value{
An object of class \code{MArrayLM} resulting from the \code{eBayes} function.
}
\details{
This function is used if the two versions of an element being compared
(often in allelic studies) are measured in the same sample. When the
\code{MPRASet} object is created, the two versions of all elements are
split into separate columns although they come from the same sample.
}
\references{
  Smyth, Gordon K., Jo\"{e}lle Michaud, and Hamish S. Scott.
  \emph{Use of within-Array Replicate Spots for Assessing Differential
    Expression in Microarray Experiments.} Bioinformatics 2005,
  21(9):2067-75. \doi{10.1093/bioinformatics/bti270}. 
}
\examples{
data(mpraSetExample)
design <- data.frame(intcpt = 1,
                     episomal = grepl("MT", colnames(mpraSetExample)))
## Create fake blocks. This creates pairings between columns
##        1 and 4, 2 and 5, 3 and 6
blocks <- rep(1:3, 2)
## Technically, this is an incorrect analysis because the columns
## are not paired
## This is just for illustration
mpralm_fit <- fit_corr(object = mpraSetExample, aggregate = "mean",
                       design = design, normalize = TRUE, block = blocks,
                       plot = FALSE)
}
