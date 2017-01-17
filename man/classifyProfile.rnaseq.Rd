\name{classifyProfile.rnaseq}
\alias{classifyProfile.rnaseq}
\title{
Expression profile classification
}
\description{
Function to classify RNA-seq gene expression profiles
}
\usage{
classifyProfile.rnaseq(ref_matrix, query_mat, gene.ids.type="ensembl", fun1 = median, write2File=FALSE, out.dir=getwd())
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
  \item{write2File}{
A logical value. If TRUE the classification results will be written to a file.
}
  \item{out.dir}{
Path to the directory, in which to write the results. Default is the actual working directory.
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




\examples{

library(sampleClassifierData)
data("se_rnaseq_refmat")
rnaseq_refmat <- assay(se_rnaseq_refmat)
data("se_rnaseq_testmat")
rnaseq_testmat <- assay(se_rnaseq_testmat)
res2.list <- classifyProfile.rnaseq(ref_matrix=rnaseq_refmat, query_mat=rnaseq_testmat, 
gene.ids.type="ensembl",write2File=FALSE)
res2.list


}
\keyword{ RNA-seq data }
\keyword{ classification }
