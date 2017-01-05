#Demo script to assign score to each miRNA##
#Debajyoti Sinha#########01 November 2016###
############################################
#Define Working Directory
setwd("ParSel-master/sc")
dir.create(file.path(getwd(), "../plot/"), showWarnings = FALSE)
dir.create(file.path(getwd(), "../output/"), showWarnings = FALSE)

#Define Cores to Use
usecores = 4


#library(RcppCNPy)
require(glmnet)
require(ROCR)
require(caret)
require(survival)
require(ggplot2)
require(reshape2)
#Call prepare script in R:
# 1. Save prepared data
# 2. Save test and train data

source("prepare.R")

#Run Clusreting R script
# input training
# Save network topology adjacency matrix
clustering<-function(){
  cat("Creating clusters...\n")
  t1 <- Sys.time()
  #source("clusters.R") #Uncomment if You choose WGCNA
  source("clusters2.R")
  t2 <- Sys.time()
  cat("Completed creating clusters.\n")
  return(cat("Clustering Time:", difftime(t2,t1),"\n"))
}

#Run python script
scoring<-function(){
  # Output gene scores
  # Start the clock!
  t1 <- Sys.time()
  system(paste("python approx_proposed_roulette_size.py", usecores))
  # Stop the clock
  t2 <- Sys.time()
  return(cat("Feature Scoring Time:",difftime(t2,t1),"\n"))
}

write.table(file="../output/gene_names.csv",gene_names,quote=FALSE,sep='\t')
variables<-substring(read.table(file="../output/gene_names.csv",sep='\t',stringsAsFactors = F)$x,7)
gene_names = read.csv("../output/gene_names.csv",header=TRUE,quote="",sep="\t");

ntimes=10
nsamples=dim(datat)[1]
set.seed(1)
splits <- createDataPartition(datat$type, p = 0.6, times = ntimes)
effect_list<-list()
for(n in 1:ntimes){
  
  nfeat<-50#length(feat_lasso)
  
  train<-t(datat[c(splits[[n]]),])
  test<-t(datat[c(setdiff(1:nsamples,splits[[n]])),])
  
  write.csv(train,file="../output/train.csv",row.names=FALSE, quote=FALSE)
  write.csv(test,file="../output/test.csv",row.names=FALSE, quote=FALSE)
  
  cat("\nBootstrap Set ",n,"...\n")

  clustering()

  scoring()

  
  if(!file.exists(paste("../output/bigfile_",n,sep = ""))) next
  our<-read.csv(paste("../output/bigfile_",n,sep = ""),sep = ',', header=F)
  our_df<-data.frame('index'=our$V1,'values'=our$V2,'effect'=our$V3)
  #ranks<-rank(-our_df$effect)
  #rank_lists<-cbind(rank_lists,ranks)
  effect<- our_df$effect
  effect_list<-append(effect_list,list(variables[intersect(order(-effect),which(effect>0))][1:nfeat]))
}


#Accent, Dark2, Paired, Pastel1, Pastel2, Set1, Set2, Set3
most_frequent<-table(as.factor(unlist(effect_list)))*100/ntimes
names(most_frequent)<-gsub("\\.", "-",  as.character(names(most_frequent)))
bar_data<-melt(sort(most_frequent[which(most_frequent>65)],decreasing = T))
pdf("../plot/top_mirna_KIRC.pdf")
ggplot(bar_data,aes(x=Var1,y = value))+geom_bar(stat="identity", width=0.8 , fill = "grey60")+theme_classic()+
  theme(axis.text.y = element_text(size=12), axis.text.x = element_text(angle = 75,vjust = 0.5, size=15),axis.line.x=element_line(size = 0.5))+
  ylab("% frequency")+xlab("Top miRNAs")+
  ggtitle("KIRC Dataset")
dev.off()
