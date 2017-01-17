### Function to classify samples 
### the function expects 3 parameters as input:
### markers_list: the list of markers returned by the marker tool
### ref_matrix: the expression matrix used by the marker tool to predict the marker genes
### query_mat: the query expression profile/s
### chip1: chip name of the reference matrix
### chip2: chip name of the query matrix
### fun1 : function (mean or median) to be used to filter the marker genes, default is median.
### fun2: function (mean or median) to be used to to summarize the expression values of probe sets that belong to the same gene. This 
#   parameter can be ignored if the reference and query matrix are from the same chip.
### write2File: If TRUE, the classification results for each query profile will be written to a file
### out.dir : Path to a directory to write the classification results, default is the current working directory

classifyProfile <- function(ref_matrix, query_mat, chip1="hgu133plus2",chip2="hgu133a", fun1=median, fun2=mean, write2File=FALSE,out.dir=getwd()){

if(dim(ref_matrix)[1]!=dim(query_mat)[1]){
cat("The reference matrix and the query are from different platforms...\n")
cat("Collapse rows ...\n")
ref.mat <- .collapseRows(ref_matrix, chip1, FUN=fun2) 
query.mat <- .collapseRows(query_mat, chip2, FUN=fun2) 
if(dim(ref.mat)[1] > dim(query.mat)[1]){
ei.inds <- match(rownames(query.mat),rownames(ref.mat))
ref.mat2 <- ref.mat[ei.inds,]
ref_matrix <- ref.mat2
query_mat <- query.mat
}
else{
ei.inds <- match(rownames(ref.mat),rownames(query.mat))
query.mat2 <- query.mat[ei.inds,]
query_mat <- query.mat2
ref_matrix <- ref.mat
}
}

#require(MGFM)
cat("detecting marker genes...\n")
mg.list <- getMarkerGenes(ref_matrix, samples2compare="all", annotate=FALSE, chip=chip1, score.cutoff=1)
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
cat("done!\n")
predicted1 <- lapply(predicted, function(x) { x["Present.markers"] <- NULL; x })
return(predicted1)
}

