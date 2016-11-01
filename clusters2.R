list.of.packages2 <- c("fastcluster")
new.packages <- list.of.packages2[!(list.of.packages2 %in% installed.packages()[,"Package"])]
if(length(new.packages)) {
  install.packages("fastcluster")
}

# Load the hclust package
require(fastcluster)
doclustering<-function(){
  #Read in the train data set
  readData = read.csv("../output/train.csv",header=TRUE)
  readData<-data.frame(as.matrix(sapply(readData, as.numeric)))
  gene_names = read.csv("../output/gene_names.csv",header=TRUE,quote="",sep="\t");
  
  m_matrix<-data.matrix(readData[-1,])
  rownames(m_matrix)<-gene_names$names
  ind <- apply(m_matrix, 1, sd) == 0
  m_matrix <- m_matrix[!ind,]
  cor_t <- 1 - abs(cor(t(m_matrix)))
  distancet <- as.dist(cor_t)#dist(m_matrix, method="minkowski")
  hclust_complete <- hclust(distancet, method = "complete")
  #dendcomplete <- as.dendrogram(hclust_complete)
  #heatmap(m_matrix, Rowv=dendcomplete, Colv=NA)
  
  nclust = 30
  memb <- cutree(hclust_complete, k = nclust)
  
  adjlist<-replicate(nclust,c())
  
  for ( i in 1:length(memb)){
    adjlist[[memb[[i]]]]<-c(adjlist[[memb[[i]]]],i)
  }
  
  
  fileConn<-file("../output/adjacency_2.csv")
  writeLines(unlist(lapply(adjlist, paste, collapse=",")), fileConn)
  close(fileConn)
}

doclustering()

#rm(list=ls(all=TRUE))