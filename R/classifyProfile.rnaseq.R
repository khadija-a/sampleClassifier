### Function to classify RNA-seq samples 
### the function expects the following parameters as input:
### ref_matrix : Reference matrix as RNA-seq gene expression matrix with genes corresponding to rows and samples corresponding to columns.
### query_mat: The query expression profile/s
### gene.ids.type: Type of the used gene identifiers, the following gene identifiers are supported: ensembl, refseq and ucsc gene ids
### fun1 : Function (mean or median) to be used to filter the marker genes, default is median.
### write2File: If TRUE, the classification results for each query profile will be written to a file
### out.dir : Path to a directory to write the classification results, default is the current working directory

classifyProfile.rnaseq <- function(ref_matrix, query_mat, gene.ids.type="ensembl",fun1=median, write2File=FALSE,out.dir=getwd()){

if(dim(ref_matrix)[1]!=dim(query_mat)[1]){
stop("The reference and the query matrix should have the same gene ids...\n")}

#require(MGFR)
x <- ref_matrix
ref_matrix <- x[!apply(x == 0, 1, all), , drop = FALSE]
mg.list <- getMarkerGenes.rnaseq(ref_matrix, samples2compare="all", gene.ids.type=gene.ids.type,annotate=FALSE,score.cutoff=1) 
markers_list <- .getmarker_ps(mg.list,fun=fun1)
if(!is.matrix(query_mat)){
query_mat <- as.matrix(query_mat)
}
inds <- match(rownames(ref_matrix), rownames(query_mat))
sort_mat <- as.matrix(query_mat[inds,])
s.len <- dim(query_mat)[2]
if(s.len > 1){
cat(s.len, "profiles to be classified...\n")
}else{
cat(s.len, "profile to be classified...\n")
}
nhits <- length(table(colnames(ref_matrix)))
predicted <- lapply(apply(sort_mat,2,FUN=.classify_sample,markers.list=markers_list, ref.matrix=ref_matrix,hits=nhits), "[[", 1)

if(write2File){
sapply(names(predicted), 
 function (x) write.csv2(predicted[[x]], file=paste(out.dir,paste(x, "csv", sep="."),sep="/")))
 }
cat("Done!\n")
predicted1 <- lapply(predicted, function(x) { x["Present.markers"] <- NULL; x })
return(predicted1)
}