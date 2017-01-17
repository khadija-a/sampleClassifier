\name{classifyProfile.rnaseq.svm}
\alias{classifyProfile.rnaseq.svm}
\title{
Expression profile classification
}
\description{
Function to classify RNA-seq gene expression profiles using support vector machines (SVM)
}
\usage{
classifyProfile.rnaseq.svm(ref_matrix, query_mat, gene.ids.type="ensembl", fun1 = median)
}
\arguments{
  \item{ref_matrix}{
RNA-seq data matrix to be used as reference, with genes corresponding to rows
and samples corresponding to columns.
}
  \item{query_mat}{
RNA-seq query matrix to be classified, with genes corresponding to rows
and samples corresponding to columns.
}
  \item{gene.ids.type}{
Type of the used gene identifiers, the following gene identifiers are supported: ensembl, refseq and ucsc gene ids. Default is ensembl.}

  \item{fun1}{
\link{mean} or \link{median}. This will specify the number of marker genes that will be used for classification. Default is \link{median}.
}

}

\value{
A data frame with the predicted classes for each query profile.
}
\details{
  This function is based on the function \link[e1071]{svm} from the R-package
'e1071'.
}
\author{
Khadija El Amrani <khadija.el-amrani@charite.de>
}




\examples{

library(sampleClassifierData)
data("se_rnaseq_refmat")
rnaseq_refmat <- assay(se_rnaseq_refmat)
data("se_rnaseq_testmat")
rnaseq_testmat <- assay(se_rnaseq_testmat)
res2.svm.df <- classifyProfile.rnaseq.svm(ref_matrix=rnaseq_refmat, query_mat=rnaseq_testmat, 
gene.ids.type="ensembl")
res2.svm.df

}
\keyword{ SVM }
\keyword{ RNA-seq data }
\keyword{ classification }
