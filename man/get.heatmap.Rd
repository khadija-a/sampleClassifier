\name{get.heatmap}
\alias{get.heatmap}
\title{
display classification results as heatmap
}
\description{
Function to display the classification predictions as a heatmap
}
\usage{
get.heatmap(res.list)
}
\arguments{
  \item{res.list}{
the result list returned by the function  \code{\link{classifyProfile}} or \code{\link{classifyProfile.rnaseq}}
}
}
\details{
  This function is based on the function \link[ggplot2]{ggplot} from the R-package
'ggplot2'.
}
\value{
This function is used only for the side effect of creating a heatmap.
}

\author{
Khadija El Amrani <khadija.el-amrani@charite.de>}

\examples{

library(sampleClassifierData)
data("se_micro_refmat")
micro_refmat <- assay(se_micro_refmat)
data("se_micro_testmat")
micro_testmat <- assay(se_micro_testmat)
res1.list <- classifyProfile(ref_matrix=micro_refmat, query_mat=micro_testmat,
chip1="hgu133plus2",chip2="hgu133a", write2File=FALSE)
get.heatmap(res1.list)
}

\keyword{ heatmap }
\keyword{ classification }
