\name{classifyProfile}
\alias{classifyProfile}
\title{
Expression profile classification
}
\description{
Function to classify microarray gene expression profiles
}
\usage{
classifyProfile(ref_matrix, query_mat, chip1 = "hgu133plus2", chip2 = "hgu133a", fun1 = median, fun2 = mean, write2File=FALSE, out.dir=getwd())
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
  \item{write2File}{
If TRUE, the classification results for each query profile will be written to a file.
}
  \item{out.dir}{
Path to a directory to write the classification results, default is the current working directory.
}
}
\details{
Each query profile is compared to all sample types in the reference matrix and a similarity score is calculated. The similarity score is based on the number of marker genes that are shared between the query and the reference. These marker genes are given in a file if write2File is TRUE.
}
\value{
A list with top hits for each query profile, sorted according to a similarity score.
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
res1.list <- classifyProfile(ref_matrix=micro_refmat, query_mat=micro_testmat,
chip1="hgu133plus2",chip2="hgu133a", write2File=FALSE)
res1.list


}
\keyword{ microarray data }
\keyword{ classification }
