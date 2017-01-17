\name{classifyProfile.svm}
\alias{classifyProfile.svm}
\title{
Expression profile classification
}
\description{
Function to classify microarray gene expression profiles using support vector machines (SVM)
}
\usage{
classifyProfile.svm(ref_matrix, query_mat, chip1 = "hgu133plus2", chip2 = "hgu133a", fun1 = median, fun2 = mean)
}
\arguments{
  \item{ref_matrix}{
Normalized microarray data matrix to be used as reference, with probe sets corresponding to rows
and samples corresponding to columns.
}
  \item{query_mat}{
  Normalized microarray query matrix to be classified, with probe sets corresponding to rows
and samples corresponding to columns.
}
  \item{chip1}{
Chip name of the reference matrix.
}
  \item{chip2}{
Chip name of the query matrix. This parameter can be ignored if the reference and query matrix are from the same chip.
}
  \item{fun1}{
\link{mean} or \link{median}. This will specify the number of marker genes that will be used for classification. Default is \link{median}.
}
  \item{fun2}{
\link{mean} or \link{median}. This will be used to summarize the expression values of probe sets that belong to the same gene. This 
parameter can be ignored if the reference and query matrix are from the same chip. Default is \link{mean}.
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



\seealso{
see also \code{\link[MGFM]{getMarkerGenes}}.
}
\examples{

library(sampleClassifierData)
data("se_micro_refmat")
micro_refmat <- assay(se_micro_refmat)
data("se_micro_testmat")
micro_testmat <- assay(se_micro_testmat)
res1.svm.df <- classifyProfile.svm(ref_matrix=micro_refmat, query_mat=micro_testmat,
chip1="hgu133plus2",chip2="hgu133a")
res1.svm.df


}
\keyword{ microarray data }
\keyword{ classification }
\keyword{ SVM }