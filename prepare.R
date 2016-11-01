#check for installed packages
list.of.packages <- c("caret", "randomForest")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)


require(randomForest)

dataRead <- function(fileName)
{
  data <- read.table(fileName, header=T,na.strings=c("[Pending]","[Not Available]","[Not Applicable]","null","null ","NA"),  sep="\t", quote="")
  samples <- as.character(data[,1])
  data <- data[,2:(ncol(data)-1)]
  rownames(data) <- samples
  cat(paste(fileName, ": ", nrow(data), "(samples),", ncol(data), "(features)\n"))
  data <- na.roughfix(data) # impute the missing value in x, require randomForest library
  valid.cols <- which(apply(data, 2, sd)> 0) # remove the non-informative ones, with same values (e.g., 0) across all samples
  cat(paste((ncol(data)-length(valid.cols)), "invalid cols were removed.\n"))
  data <- data[,valid.cols]    
  return (data)
}

cat("Reading raw data matrix...\n")
s<-Sys.time()
mRNA_core<-dataRead("../data/KIRC_miRNA_core.txt")
d_t<-Sys.time() -s
binary_survival <- read.delim("../data/KIRC_binary_survival.txt")
match_ids<-match(binary_survival$feature,rownames(mRNA_core))
shape=dim(mRNA_core)
datat<-data.frame(cbind("type"=binary_survival$is_alive,mRNA_core[match_ids,]))
gene_names<-data.frame("names"=toupper(unlist(lapply(X = colnames(mRNA_core),FUN = function(X){unlist(strsplit(X,"[_. ]"))[2]}))))
#train.all <- read.table("../data/KIRC_train_sample_list.txt",header=F, stringsAsFactors=F)
#test.all <- read.table("../data/KIRC_test_sample_list.txt",header=F, stringsAsFactors=F)


gene_names<-colnames(mRNA_core)
surv.data <-dataRead("../data/KIRC_OS_core.txt") # survival data
mySurv <- Surv(surv.data$OS_OS, surv.data$OS_vital_status,type='right')
rownames(mySurv)<-rownames(mRNA_core)
print(paste("Total events:", length(which(mySurv[,"status"]==1))))

#rm(list=ls(all=TRUE))
