\name{mpraSetExample}
\alias{mpraSetExample}
\alias{mpraSetAggExample}
\alias{mpraSetAllelicExample}
\docType{data}
\title{
  Example data for the mpra package.
}
\description{
  Example data for the MPRA package. \code{mpraSetExample} and 
  \code{mpraSetAggExample} come from a study by Inoue et al 
  that compares episomal and lentiviral MPRA. The former contains 
  data at the barcode level and the latter contains data 
  aggregated over barcodes. \code{mpraSetAllelicExample} come from 
  a study by Tewhey et al that looks at regulatory activity of 
  allelic versions of thousands of SNPs to follow up on prior 
  eQTL results.
}
\usage{
data("mpraSetExample")
data("mpraSetAggExample")
data("mpraSetAllelicExample")
}
\format{
  An \code{MPRASet}.
}
\details{
  \code{mpraSetExample} contains barcode level information for the 
  study by Inoue et al.
  \code{mpraSetAggExample} contains count information from 
  \code{mpraSetExample} where the counts have been summed over 
  barcodes for each element.
  \code{mpraSetAllelicExample} contains count information for the 
  Tewhey et al study. The counts have been summed over barcodes 
  for each element.
}
\source{
  A script for creating the three datasets is supplied in the
  \code{scripts} folder of the package. The data are taken from the GEO
  submission associated with the paper (see references), specifically
  GSE83894 and GSE75661.
}
\references{
  Inoue, Fumitaka, Martin Kircher, Beth Martin, Gregory M. Cooper, 
  Daniela M. Witten, Michael T. McManus, Nadav Ahituv, and 
  Jay Shendure. \emph{A Systematic Comparison Reveals Substantial 
  Differences in Chromosomal versus Episomal Encoding of Enhancer 
  Activity}. Genome Research 2017, 27(1):38-52.
  \doi{10.1101/gr.212092.116}.

  Tewhey R, Kotliar D, Park DS, Liu B, Winnicki S, Reilly SK, Andersen KG, Mikkelsen TS, Lander ES, Schaffner SF, Sabeti PC. \emph{Direct Identification of Hundreds of Expression-Modulating Variants using a Multiplexed Reporter Assay}. Cell 2016, 165:1519-1529. \doi{10.1016/j.cell.2016.04.027}.
}
\examples{
data(mpraSetAggExample)
}
\keyword{datasets}
