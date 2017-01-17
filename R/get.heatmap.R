get.heatmap <- function(res.list){
#require(ggplot2)
Samples <- Reference <- Score <- NULL
nsamples <- length(res.list)
nhits <- dim(res.list[[1]])[1]
nmat <- matrix(0, ncol=nhits,nrow=nsamples)
colnames(nmat) <- res.list[[1]][,1][order(res.list[[1]][,1])]
rownames(nmat) <- names(res.list)
for(i in 1:nsamples){
inds <- match(colnames(nmat),as.character(res.list[[i]][,1]))
nmat[i,] <- as.numeric(as.character(res.list[[i]][inds,2]))
}
mydf <- data.frame()
for(i in 1:nsamples){
temp.df <- as.data.frame(cbind(rep(rownames(nmat)[i],nhits),colnames(nmat),nmat[i,]), stringsAsFactors=FALSE)
mydf <- rbind(mydf,temp.df)
}
mydf[,3] <- as.numeric(mydf[,3])
colnames(mydf) <- c("Samples","Reference","Score")
base_size <-10
    (p <- ggplot(mydf, aes(Reference, Samples)) +
    geom_tile(aes(fill = Score), colour = "black") +
    scale_fill_gradient(low = "black",high = "green",limits=c(0,1),breaks=seq(0, 1, by=0.2)))
    
     p  + labs(x = "Reference",
     y = "Query samples")+ theme( axis.text.x = element_text(size = base_size *
        0.8, angle = 330, hjust = 0, colour = "grey50"))
}


