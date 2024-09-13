\name{mpra-package}
\alias{mpra-package}
\alias{mpra}
\docType{package}
\title{
\packageTitle{mpra}
}
\description{
\packageDescription{mpra}
}
\details{
  This package provides tools for the analysis of MPRA data. The primary
  purpose is to enable powerful differential analysis of activity
  measures, but the package can also be used to generate precision
  weights useful in regression analyses of activity scores on sequence
  features. The main workhorse is the \code{mpralm} function which draws
  on the previously proposed \code{voom} framework for RNA-seq analysis
  in the \code{limma} package.
}
\author{
\packageAuthor{mpra}

Maintainer: \packageMaintainer{mpra}
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
\keyword{ package }
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
