\name{normalize_counts}
\alias{normalize_counts}
\title{Total count normalization of DNA and RNA counts}
\description{
Total count normalization of DNA and RNA counts.
}
\usage{
normalize_counts(object, block = NULL)
}
\arguments{
  \item{object}{An object of class \code{MPRASet}.}
  \item{block}{
A vector giving the sample designations of the columns of \code{object}. The default, \code{NULL}, indicates that all columns are separate samples.
}
}
\details{
\code{block} is a vector that is used when the columns of the \code{MPRAset} object are paired. This often is the case when comparing allelic versions of an element. In this case, the first $s$ columns of \code{object} give the counts for the reference allele in $s$ samples. The second $s$ columns give the counts for the alternative allele measured in the same $s$ samples. With 3 samples, \code{block} would be \code{c(1,2,3,1,2,3)}. All columns are scaled to have 10 million counts.
}
\value{
An object of class \code{MPRASet} with the total count-normalized DNA and RNA counts.
}
\examples{
data(mpraSetAggExample)
mpraSetAggExample <- normalize_counts(mpraSetAggExample)
}
