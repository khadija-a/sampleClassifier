\name{sampleClassifier-package}
\alias{sampleClassifier-package}
\alias{sampleClassifier}
\docType{package}
\title{
Sample Classifier
}
\description{
The package is designed to classify samples from microarray and RNA-seq gene expression datasets.
}
\details{
\tabular{ll}{
Package: \tab sampleClassifier\cr
Type: \tab Package\cr
Version: \tab 1.0.0\cr
License: \tab GPL-3 \cr
}

}
\author{
Khadija El Amrani
Maintainer: Khadija El Amrani <khadija.el-amrani@charite.de>}

\keyword{ package }

\examples{
\dontrun{
library(sampleClassifierData)
data("se_micro_refmat")
micro_refmat <- assay(se_micro_refmat)
data("se_micro_testmat")
micro_testmat <- assay(se_micro_testmat)
res1.list <- classifyProfile(ref_matrix=micro_refmat, query_mat=micro_testmat,
chip1="hgu133plus2",chip2="hgu133a", write2File=FALSE)
res1.list
}
}
