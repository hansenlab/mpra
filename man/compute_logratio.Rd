\name{compute_logratio}
\alias{compute_logratio}
\title{
  Compute activity measure (log-ratio) for each element.
}
\description{
Compute the log ratio of RNA counts to DNA counts using different methods. 
For \code{"mean"}, uses the average of barcode-specific log ratios. 
For \code{"sum"}, sums RNA and DNA counts over barcodes before forming the log ratio.
}
\usage{
compute_logratio(object, aggregate = c("mean", "sum", "none"))
}
\arguments{
  \item{object}{An object of class \code{MPRASet}.}
  \item{aggregate}{Aggregation method over barcodes: \code{"mean"} to
    use the average of barcode-specific log ratios, \code{"sum"} to use
    the log ratio of summed RNA and DNA counts, \code{"none"} to perform
    no aggregation (counts have already been summarized over barcodes).}
}
\value{
A \code{matrix} with the same dimension as \code{object}, containing element- 
and sample-specific log ratios.
}

\examples{
  data(mpraSetAggExample)
  logr <- compute_logratio(mpraSetAggExample, aggregate = "sum")
}
