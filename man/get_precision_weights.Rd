\name{get_precision_weights}
\alias{get_precision_weights}
\title{Get precision weights from the mean-variance relationship.}
\description{
Estimates the variability of the supplied log-ratios across samples as a function of copy number (DNA count levels).
}
\usage{
get_precision_weights(logr, design, log_dna, span = 0.4, plot = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{logr}{
Matrix of outcome measures: log2 ratio of RNA counts to DNA counts.
}
  \item{design}{
Design matrix specifying comparisons of interest.
}
  \item{log_dna}{
Matrix of log2 aggregated DNA counts of the same dimension as \code{logr}.
}
  \item{span}{
The smoothing span for \code{lowess} in estimating the mean-variance relationship. Default: 0.4.
}
  \item{plot}{
If \code{TRUE}, plot the mean-variance relationship.
}
  \item{\dots}{
Further arguments to be passed to \code{lmFit} for obtaining residual standard deviations used in estimating the mean-variance relationship.
}
}
\details{
Residual standard deviations are computed using the supplied outcomes and design matrix. The square root of the the residual standard deviations are modeled as a function of the average log2 aggregated DNA counts to estimate the mean-variance relationship.
}
\value{
A matrix of precision weights of the same dimension as \code{logr} and \code{log_dna}.
}
\references{
  Law, Charity W., Yunshun Chen, Wei Shi, and Gordon K. Smyth. 2014. 
  \emph{Voom: Precision Weights Unlock Linear Model Analysis Tools for RNA-Seq Read Counts}. 
  Genome Biology 15 (2): R29. doi:10.1186/gb-2014-15-2-r29.
}
\examples{
data(mpraSetAggExample)
design <- data.frame(intcpt = 1,
                     episomal = grepl("MT", colnames(mpraSetAggExample)))
logr <- compute_logratio(mpraSetAggExample, aggregate = "none")
log_dna <- log2(getDNA(mpraSetAggExample, aggregate = FALSE) + 1)
w <- get_precision_weights(logr = logr, design = design,
                           log_dna = log_dna, plot = FALSE)
}
