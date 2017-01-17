### Function to classify a single sample 
### this an intern function which will be called from the function classify_profile_opt
### the function expects 3 parameters as input:
### markers.list: the list of markers returned by the marker tool
### ref.matrix: the expression matrix used by the marker tool to predict the marker genes
### query_profile: the query expression profile
### number of hits to return, default is 5
.get.genes <- function(IDs=c(''),chip=NULL){
#require(annotate)
 if (!require(paste(chip, "db", sep="."),character.only = TRUE))
    {
      install.packages(paste(chip, "db", sep="."),dependencies=TRUE)
        if(!require(paste(chip, "db", sep="."),character.only = TRUE)) stop("Package not found")
    }
    gs.list <- lookUp(IDs,paste(chip, "db", sep="."), "SYMBOL")
    egids.list <-lookUp(IDs,paste(chip, "db", sep="."), "ENTREZID")
    un.gs <- unlist(gs.list)
    ps.inds <- match(IDs,names(un.gs))
    un.egids <- unlist(egids.list)
    ps.inds1 <- match(IDs,names(un.egids))
    annot.df <- as.data.frame(cbind(IDs,unname(un.gs[ps.inds]), unname(un.egids[ps.inds1])))
    colnames(annot.df) <- c("ID", "Gene Symbol","ENTREZ_GENE_ID")
	return(annot.df)
}
.classify_sample <- function(markers.list, ref.matrix, query_profile,hits=5){
tissues <- c()
scores <- c()
markers.len <- c()
new.mat <- cbind(ref.matrix, query_profile)
colnames(new.mat) <- c(colnames(ref.matrix),"query")
query.ind <- length(colnames(new.mat))
markers.vec <- c()
for(i in 1:length(markers.list)){
tissue.name <- unlist(strsplit(names(markers.list)[i],split="_markers"))
inds <- grep(tissue.name,colnames(ref.matrix),fixed=TRUE)
num.t <- length(inds)+1
ord.mat <- apply(new.mat[markers.list[[i]],,drop=FALSE], 1,order,decreasing=TRUE)
ind.m <- which(ord.mat[1:num.t,]==query.ind,arr.ind=TRUE)
pres.markers <- colnames(ord.mat)[ind.m[,2]]
tissue.score <- length(pres.markers)
tissues <- c(tissues, tissue.name)
scores <- c(scores, tissue.score)
markers.vec <- c(markers.vec, paste(pres.markers,collapse=", "))
markers.len <- c(markers.len, length(markers.list[[i]]))
} # end for loop
ratio.vec <- scores/markers.len
max.score <- max(ratio.vec)
if(max.score > 0){
ratio.out <- paste(scores,markers.len,sep=" / ")
names(ratio.out) <- tissues
names(ratio.vec) <- tissues
names(markers.vec) <- tissues
sorted.ratio <- sort(ratio.vec, decreasing=TRUE)
mhits <- min(hits,length(sorted.ratio>0))
names.inds <- match(names(sorted.ratio), names(ratio.out))
sorted.ratio.out <- ratio.out[names.inds]
sorted.markers <- markers.vec[names.inds]
out.df <-as.data.frame(cbind(names(sorted.ratio)[1:mhits],round(as.numeric(sorted.ratio)[1:mhits],digits=3), as.character(sorted.ratio.out)[1:mhits],as.character(sorted.markers)[1:mhits]),stringsAsFactors=FALSE)
colnames(out.df) <- c("Hits","Score","Ratio","Present.markers")
return(list(out.df))
}else{
return(list(""))
}
}



.getmarker_ps <- function(markers_list,fun=median){
res.list <- list()
names.vec <- names(markers_list)
res.names <- c()
lim <- 0
len.mar <- 0
llen <- c()
for(i in 1:length(markers_list)){
if(length(markers_list[[i]])>0){
llen <- c(llen, length(markers_list[[i]]))
}
}
len.mar <- round(fun(llen))
for(i in 1:length(markers_list)){
if(length(markers_list[[i]])>0){
lim <- min(len.mar, length(markers_list[[i]]))
ulist <- unlist(strsplit(markers_list[[i]][1:lim], split=" : "))
mar_ps <- ulist[seq(1,length(ulist),by=2)]
res.list[[length(res.list)+1]] <- mar_ps
res.names <- c(res.names,names.vec[i])
}
}
names(res.list) <- res.names
return(res.list)
}
.collapseRows <- function(data.mat, chip, FUN=mean){
#assign("ann.df", get(paste(chip,"df",sep=".")))
ann.df <- .get.genes(rownames(data.mat),chip=chip)
colnames.unique <- TRUE
if(length(colnames(data.mat))!=length(unique(colnames(data.mat)))){
colnames(data.mat) <- make.unique(colnames(data.mat),sep="|")
colnames.unique <- FALSE
}
entrez.ind <- grep("entrez",colnames(ann.df),ignore.case=TRUE)
amb.inds <-grep(",", ann.df[,entrez.ind])
uann.df <- ann.df
if(length(amb.inds) > 0){
uann.df <- ann.df[-amb.inds,] 
}
miss.inds <- which(is.na(uann.df[,entrez.ind]))
if(length(miss.inds)>0){
uann.df <- uann.df[-miss.inds,] 
}
ps <- as.character(uann.df[,1])
ps.inds <- match(ps,rownames(data.mat))
data.df <- cbind.data.frame(data.mat[ps.inds,],Entrez.Gene=as.factor(uann.df[,entrez.ind]),stringsAsFactors = FALSE)
colnames(data.df) <- c(colnames(data.mat), "Entrez.Gene")
res.df <- aggregate(. ~ Entrez.Gene, data = data.df, FUN)
res.df2 <- res.df[,2:dim(res.df)[2]]
rownames(res.df2) <- as.character(res.df[,1])
if(!colnames.unique){
split.list <- strsplit(colnames(res.df2), split = "[|]")
colnames(res.df2) <- sapply(split.list,function(x) x[1])
}
res.mat <- as.matrix(res.df2)
return(res.mat)
}

#### function to get matrix with marker probe sets only
#### the function expects the following parameters:
#### orig_mat: expression matrix from which to select the probe sets. The probe sets correspond to the rows and the samples to the columns
#### markers_list: list with markers of each tissue in an entry. Please note that the list of markers should be an output of the marker tool (called with annotate=TRUE)
#### fun: one of:
#### "mean" to get a submatrix with a number of probe sets that corresponds to min(mean(# of all markers), length(markers of a tissue))
#### "median" to get a submatrix with a number of probe sets that corresponds to min(median(# of all markers), length(markers of a tissue))
#### may be this is still not clear what is meant by fun="mean" or fun="median", I will try to explain it using an example
#### let us say we have a marker list with 3 entries: the first entry contains 5 markers, the 2. entry contains 10 and the 3. entry contains 15 marker probe sets
#### if fun="mean" the result matrix will contain 5 marker probe sets from the first entry (because min(5, mean(5,10,15))=5), 10 markers from the second entry ((because min(10, mean(5,10,15))=10))
##### and 10 markers from the third entry ((because min(15, mean(5,10,15))=10))
##### the same is true for median
.get_mpsmat <- function(orig_mat, markers_list, fun=median){
ps <- c()
llen <- c()
for(i in 1:length(markers_list)){
if(length(markers_list[[i]])>0){
llen <- c(llen, length(markers_list[[i]]))
}
}
mar.len <- round(fun(llen))
for(i in 1:length(markers_list)){
if(length(markers_list[[i]])>0){
lim <- min(mar.len, length(markers_list[[i]]))
ulist <- unlist(strsplit(markers_list[[i]][1:lim],split=" : "))
ps_temp <- ulist[seq(1,length(ulist),by=2)]
ps <- c(ps, ps_temp)
}

}
ps.inds <- match(ps, rownames(orig_mat))
res_mat <- orig_mat[ps.inds,]

return(res_mat)
}