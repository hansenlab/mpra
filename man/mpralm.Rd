\name{mpralm}
\alias{mpralm}
\title{Linear models for differential analysis of MPRA data}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mpralm(object, aggregate = c("mean", "sum", "none"), normalize = TRUE, design, block = NULL, model_type = c("indep_groups", "corr_groups"), ...)
}
\arguments{
  \item{object}{An object of class \code{MPRASet}.}
  \item{aggregate}{
Aggregation method over barcodes: \code{"mean"} to use the average of barcode-specific log ratios, \code{"sum"} to use the log ratio of summed RNA and DNA counts, \code{"none"} to perform no aggregation (counts have already been summarized over barcodes).
}
  \item{normalize}{
If \code{TRUE}, perform total count normalization before model fitting.
}
  \item{design}{
Design matrix specifying comparisons of interest.
}
  \item{block}{
A vector giving the sample designations of the columns of \code{object}. The default, \code{NULL}, indicates that all columns are separate samples.
}
  \item{model_type}{
Indicates whether an unpaired model fit (\code{"indep_groups"}) or a paired mixed-model fit ((\code{"corr_groups"})) should be used.
}
  \item{\dots}{
Further arguments to be passed to \code{lmFit} for obtaining residual standard deviations used in estimating the mean-variance relationship.
}
}
\details{
Using \code{method_type = "corr_groups"} use the \code{duplicateCorrelation} function from the \code{limma} package to estimate the intra-replicate correlation of log-ratio values.
}
\value{
An object of class \code{MArrayLM} resulting from the \code{eBayes} function.
}
\references{
  Law, Charity W., Yunshun Chen, Wei Shi, and Gordon K. Smyth. 2014. 
  \emph{Voom: Precision Weights Unlock Linear Model Analysis Tools for RNA-Seq Read Counts}. 
  Genome Biology 15 (2): R29. doi:10.1186/gb-2014-15-2-r29.

  Smyth, Gordon K., Jo\"{e}lle Michaud, and Hamish S. Scott. 2005.
  \emph{Use of within-Array Replicate Spots for Assessing Differential Expression in Microarray Experiments.} Bioinformatics  21 (9): 2067-75. doi:10.1093/bioinformatics/bti270.
}
\examples{
  data(mpraSetAggExample)
  design <- data.frame(intcpt = 1, episomal = grepl("MT", colnames(getDNA(mpraSetAggExample))))
  mpralm_fit <- mpralm(object = mpraSetAggExample, aggregate = "none", normalize = TRUE, design = design, model_type = "indep_groups")
  toptab <- topTable(fit, coef = 2, number = Inf)
  head(toptab)
}
