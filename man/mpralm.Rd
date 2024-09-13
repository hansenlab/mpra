\name{mpralm}
\alias{mpralm}
\title{Linear models for differential analysis of MPRA data}
\description{
Fits weighted linear models to test for differential activity in MPRA data.
}
\usage{
mpralm(object, design, aggregate = c("mean", "sum", "none"),
       normalize = TRUE, normalizeSize = 10e6,
       block = NULL, model_type = c("indep_groups", "corr_groups"),
       plot = TRUE, endomorphic = FALSE, ...)
}
\arguments{
  \item{object}{An object of class \code{MPRASet}.}
  \item{design}{Design matrix specifying comparisons of interest. The 
    number of rows of this matrix should equal the number of columns 
    in \code{object}. The number of columns in this design matrix has 
    no constraints and should correspond to the experimental design.}
  \item{aggregate}{Aggregation method over barcodes: \code{"mean"} to
    use the average of barcode-specific log ratios, \code{"sum"} to use
    the log ratio of summed RNA and DNA counts, \code{"none"} to perform
    no aggregation (counts have already been summarized over barcodes).}
  \item{normalize}{If \code{TRUE}, perform total count normalization
    before model fitting.}
  \item{normalizeSize}{If normalizing, the target library size (default
    is 10e6).}
  \item{block}{A vector giving the sample designations of the columns of
    \code{object}. The default, \code{NULL}, indicates that all columns
    are separate samples.} 
  \item{model_type}{Indicates whether an unpaired model fit
    (\code{"indep_groups"}) or a paired mixed-model fit
    ((\code{"corr_groups"})) should be used.}
  \item{plot}{If \code{TRUE}, plot the mean-variance relationship.}
  \item{endomorphic}{If \code{TRUE}, return the same class as the input,
    i.e. an object of class \code{MPRASet}.}
  \item{\dots}{Further arguments to be passed to \code{lmFit} for
    obtaining residual standard deviations used in estimating the
    mean-variance relationship.} 
}
\details{
Using \code{method_type = "corr_groups"} use the
\code{duplicateCorrelation} function from the \code{limma} package to
estimate the intra-replicate correlation of log-ratio values.
}
\value{
An object of class \code{MArrayLM} resulting from the \code{eBayes}
function.

If \code{endomorphic = TRUE}, then an \code{MPRASet} is returned,
with the output of \code{topTable} added to the \code{rowData},
and the \code{MArrayLM} results added as an attribute
\code{"MArrayLM"}.
}
\references{
  Myint, Leslie, Dimitrios G. Avramopoulos, Loyal A. Goff, and Kasper
  D. Hansen.
  \emph{Linear models enable powerful differential activity analysis in
    massively parallel reporter assays}. 
  BMC Genomics 2019, 209. \doi{10.1186/s12864-019-5556-x}.

  Law, Charity W., Yunshun Chen, Wei Shi, and Gordon K. Smyth.
  \emph{Voom: Precision Weights Unlock Linear Model Analysis Tools for RNA-Seq Read Counts}. 
  Genome Biology 2014, 15:R29. \doi{10.1186/gb-2014-15-2-r29}.

  Smyth, Gordon K., Joelle Michaud, and Hamish S. Scott. 
  \emph{Use of within-Array Replicate Spots for Assessing Differential
    Expression in Microarray Experiments.} 
  Bioinformatics 2005, 21 (9): 2067-75. \doi{10.1093/bioinformatics/bti270}.
}
\examples{
data(mpraSetAggExample)
design <- data.frame(intcpt = 1,
                     episomal = grepl("MT", colnames(mpraSetAggExample)))
mpralm_fit <- mpralm(object = mpraSetAggExample, design = design,
                     aggregate = "none", normalize = TRUE, 
                     model_type = "indep_groups", plot = FALSE)
toptab <- topTable(mpralm_fit, coef = 2, number = Inf)
head(toptab)
}
