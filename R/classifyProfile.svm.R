classifyProfile.svm <- function(ref_matrix, query_mat, chip1="hgu133plus2",chip2="hgu133a", fun1=median, fun2=mean){
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

