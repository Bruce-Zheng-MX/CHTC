#library("seqinr")
#library("ngramr")#
#library("phylogram")
#library('rmatio')
#library(R.matlab)
library(igraph)
#library(ecp)      
library(energy)
#library(networkD3)
#library(GreedySBTM)
#require(ade4) 
#require("gSeg")
#setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("objTest_fctns.R") 
source("binary_seg.R")

d<-get(load("reality_mining_1392.RData"))
G_list<-list()
for (i in 1:(dim(d)[3]/6)){
  G_list[[i]]<-rowSums(d[,,(i*6-5):(i*6)],dims=2)
}

G_graph<-c()
for (i in 1:(dim(d)[3]/6)){
  m<-c(as.matrix(as.matrix(G_list[[i]], "adjacency")))
  G_graph<-rbind(G_graph,m)
}
G_graph[G_graph >=1] <- 1
#sum(G_graph[1,])
#dim(G_graph)#232 
Data<-G_graph
n<-nrow(Data)

iter<-as.numeric(commandArgs(TRUE)[1])
type<-as.numeric(commandArgs(TRUE)[2])
L<-10;minLen<-10
c<-0.2;num_permut<-1000
seed_I<-get(load("seed_I_1_85_minLen10.Rdata"))
rand_I<-get(load("seed_I_105_232_minLen10.Rdata"))
#iterate every interval
max_stat_index<-c();p_val_interval<-c();critical_value_interval<-c()


#for (i in 1:nrow(rand_I)){
#for (i in 1:nrow(seed_I)){
  #construct distmat for each interval
  #data<-Data[rand_I[i,1]:rand_I[i,2],]
  i<-iter
  if (type==1){
  data<-Data[seed_I[i,1]:seed_I[i,2],]
}else {data<-Data[rand_I[i,1]:rand_I[i,2],]}
  #data<-Data[seed_I[i,1]:seed_I[i,2],]
  n<-nrow(data)

  max_stat<-c();
  
  for (j in 0:num_permut){
    cat(i,'th interval',j,"permutation","\n")
    if (j!=0){
      data_permut<-data[sample(nrow(data)),]
    }else{data_permut<-data}
  distmat<-as.matrix(dist(data_permut , method = 'euclidean' ))
  
  test<-c()
  for (cp in seq(n*c,n*(1-c),1)){
    seg1<-data_permut[1:cp,];seg2<-data_permut[(cp+1):n,]
    
    if (cp<=1){
      testStat <- getT( distmat = distmat, indices = 1:n, n = 1, m = nrow(seg2), cut_off = cut_off )}
    
    else if (cp>=n-1){
      testStat <- getT( distmat = distmat, indices = 1:n, n = nrow(seg1), m = 1, cut_off = cut_off )}
    
    else {testStat <- getT( distmat = distmat, indices = 1:n, n = nrow(seg1), m = nrow(seg2), cut_off = cut_off )}
    test<-cbind(test,testStat)
  }


  max_stat<-cbind(max_stat,apply(test,1,max))
  if (j==0){
    max_stat_index<-rbind(max_stat_index,apply(test,1,which.max)+as.integer(n*c)-1)}
  
  #sum(max_stat[,1]>max_stat[,2:ncol(max_stat)])
  #dec<-as.integer(max_stat[,1]>critical_value)
  }
  critical_value<-apply(max_stat[,2:ncol(max_stat)],1,quantile,probs=0.95)
  p_val<-rowSums(max_stat[,1]<max_stat[,2:ncol(max_stat)])/num_permut
  
  #p_val_interval<-rbind(p_val_interval,p_val)
  #critical_value_interval<-rbind(critical_value_interval,critical_value)
 # }

result<-list()
result[[1]]<-p_val
result[[2]]<-max_stat_index
result[[3]]<-critical_value
result[[4]]<-max_stat
path<-paste("iter_",iter,'_type_',type,'.Rdata',sep="")
save(result,file=path)
