list.of.packages2 <- c("WGCNA")
new.packages <- list.of.packages2[!(list.of.packages2 %in% installed.packages()[,"Package"])]
if(length(new.packages)) {
  source("http://bioconductor.org/biocLite.R") 
  biocLite("WGCNA")
}

# Load the WGCNA package
require(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
enableWGCNAThreads()
#Read in the train data set
readData = read.csv("../output/train.csv",header=TRUE)
readData<-data.frame(as.matrix(sapply(readData, as.numeric)))
gene_names = read.csv("../output/gene_names.csv",header=TRUE,quote="",sep="\t");
# Take a quick look at what is in the data set:
dim(readData);
names(readData);
nSets=1
multiExpr = vector(mode = "list", length = nSets)

multiExpr[[1]] = list(data = as.data.frame(data.matrix(t(readData[-1,which(readData[1,] == 1)]))))
names(multiExpr[[1]]$data) = gene_names$x;
rownames(multiExpr[[1]]$data) = paste(readData[1,which(readData[1,] == '1')],names(readData[-1,which(readData[1,] == '1')]),sep="_");
exprSize = checkSets(multiExpr)


# Choose a set of soft-thresholding powers
powers = c(seq(4,10,by=1), seq(12,20, by=2));
# Initialize a list to hold the results of scale-free analysis
powerTables = vector(mode = "list", length = nSets);
# Call the network topology analysis function for each set in turn
for (set in 1:nSets)
  powerTables[[set]] = list(data = pickSoftThreshold(multiExpr[[set]]$data, powerVector=powers,
                                                     verbose = 2)[[2]]);
collectGarbage();

net = blockwiseConsensusModules(
  multiExpr, power = 6, minModuleSize = 30, deepSplit = 2,
  pamRespectsDendro = FALSE,
  mergeCutHeight = 0.25, numericLabels = TRUE,
  minKMEtoStay = 0, maxBlockSize = 10000,
  saveTOMs = TRUE, verbose = 2)

adjlist<-replicate(length(unique(net$color)),list())

for ( i in 1:length(net$color)){
  adjlist[[net$color[i]+1]]<-c(adjlist[[net$color[i]+1]],i)
}
fileConn<-file("../output/adjacency.csv")
writeLines(unlist(lapply(adjlist, paste, collapse=",")), fileConn)
close(fileConn)

#rm(list=ls(all=TRUE))