classifyProfile.rnaseq.svm <- function(ref_matrix, query_mat, gene.ids.type="ensembl", fun1=median){

if(dim(ref_matrix)[1]!=dim(query_mat)[1]){
stop("The reference and the query matrix should have the same gene ids...\n")}

#require(MGFR)
if(mode(ref_matrix)!="numeric"){
mode(ref_matrix) <- "numeric"
}
if(mode(query_mat)!="numeric"){
mode(query_mat) <- "numeric"
}
x <- ref_matrix
ref_matrix <- x[!apply(x == 0, 1, all), , drop = FALSE]
mg.list <- getMarkerGenes.rnaseq(ref_matrix, samples2compare="all", gene.ids.type=gene.ids.type,annotate=FALSE, score.cutoff=1)
markers_list <- .getmarker_ps(mg.list,fun=fun1)
#require(e1071)
ref.mat.svm <- ref_matrix
colnames(ref.mat.svm) <- make.unique(colnames(ref_matrix))
mpsmat <- .get_mpsmat(ref.mat.svm, mg.list, fun=fun1)
cat("building an SVM model...\n")
svm.model <- svm(x=t(mpsmat), y=colnames(ref_matrix), type="C-classification", kernel="linear")
svm.pred.class <- c()
query.mat.svm <- .get_mpsmat(query_mat, mg.list, fun=fun1)
query.mat <-as.matrix(query.mat.svm)
row.inds <- match(rownames(mpsmat),rownames(query.mat))
query.mat <- query.mat[row.inds,]
#print(table(rownames(query.mat)==rownames(mpsmat)))
#### predict classes of test matrix
s.len <- dim(query_mat)[2]
if(s.len > 1){
cat(s.len, "profiles to be classified...\n")
}else{
cat(s.len, "profile to be classified...\n")
}
predicted.svm <- predict(svm.model, t(query.mat))
predicted.svm.df <- as.data.frame(cbind(query_name=names(predicted.svm), predicted_class=as.character(predicted.svm)))
cat("done!\n")
return(predicted.svm.df)
}