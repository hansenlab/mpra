\name{compute_logratio}
\alias{compute_logratio}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
compute_logratio(object, aggregate = c("mean", "sum", "none"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An object of class \code{MPRASet}.}
  \item{aggregate}{
Aggregation method over barcodes: \code{"mean"} to use the average of barcode-specific log ratios, \code{"sum"} to use the log ratio of summed RNA and DNA counts, \code{"none"} to perform no aggregation (counts have already been summarized over barcodes).
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (object, aggregate = c("mean", "sum", "none")) 
{
    aggregate <- match.arg(aggregate)
    if (aggregate == "sum") {
        dna <- getDNA(object, aggregate = TRUE)
        rna <- getRNA(object, aggregate = TRUE)
        logr <- log2(rna + 1) - log2(dna + 1)
    }
    else if (aggregate == "mean") {
        dna <- getDNA(object, aggregate = FALSE)
        rna <- getRNA(object, aggregate = FALSE)
        eid <- getEid(object)
        logr <- log2(rna + 1) - log2(dna + 1)
        by_out <- by(logr, eid, colMeans, na.rm = TRUE)
        logr <- do.call("rbind", by_out)
        rownames(logr) <- names(by_out)
    }
    else if (aggregate == "none") {
        dna <- getDNA(object, aggregate = FALSE)
        rna <- getRNA(object, aggregate = FALSE)
        logr <- log2(rna + 1) - log2(dna + 1)
    }
    return(logr)
  }
}
