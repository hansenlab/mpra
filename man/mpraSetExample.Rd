\name{mpraSetExample}
\alias{mpraSetExample}
\alias{mpraSetAggExample}
\docType{data}
\title{
  Example data for the mpra package.
}
\description{
  Example data for the MPRA package. The package contains data both at
  the barcode level and aggregated over barcodes. The data comes from a
  study by Inoue et al that compares episomal and lentiviral MPRA.
}
\usage{
data("mpraSetExample")
data("mpraSetAggExample")
}
\format{
  An \code{MPRASet}.
}
\details{
  \code{mpraSetExample} contains barcode level information. 
  \code{mpraSetAggExample} contains count information where the counts 
  have been summed over barcodes for each element.
}
\source{
  FIMXE
}
\references{
  Inoue, Fumitaka, Martin Kircher, Beth Martin, Gregory M. Cooper, 
  Daniela M. Witten, Michael T. McManus, Nadav Ahituv, and 
  Jay Shendure. 2017. \emph{A Systematic Comparison Reveals Substantial 
  Differences in Chromosomal versus Episomal Encoding of Enhancer 
  Activity}. Genome Research 27 (1): 38–52. doi:10.1101/gr.212092.116.
}
\examples{
data(mpraSetAggExample)
}
\keyword{datasets}
