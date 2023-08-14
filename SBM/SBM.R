
#setwd(dirname(rstudioapi::getSourceEditorContext()$path))


source("objTest_fctns.R") 
source("binary_seg.R") 
#require(sbm)
require(igraph)
#require(ecp)

iter<-as.numeric(commandArgs(TRUE)[1])
L<-10;minLen<-10
c<-0.1;num_permut<-1000

seed_I<-seeded.intervals(400,decay = sqrt(2))
seed_I<-seed_I[(seed_I[,2]-seed_I[,1])>minLen,]
seed_I<-seed_I[seed_I[,2]<100|seed_I[,1]>100,]
seed_I<-seed_I[seed_I[,2]<200|seed_I[,1]>200,]
nrow(seed_I)


max_stat_index<-c();p_val_interval<-c();critical_value_interval<-c()
type<-1
i=iter

Data<-get(load("SBM_data_nodes200.Rdata"))
if (type==1){
  data<-Data[seed_I[i,1]:seed_I[i,2],]
}else {data<-Data[rand_I[i,1]:rand_I[i,2],]}
n<-nrow(data)
max_stat<-c();

for (j in 0:100){
  cat(i,'th interval',j,"permutation","\n")
  if (j!=0){
    data_permut<-data[sample(nrow(data)),]
  }else{data_permut<-data}
  distmat<-as.matrix(dist(data_permut , method = 'euclidean' ))
  
  test<-c()
  for (cp in seq(max(2,n*c),min(n*(1-c),n-2),1)){
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
  
}
critical_value<-apply(max_stat[,2:ncol(max_stat)],1,quantile,probs=0.95)
p_val<-rowSums(max_stat[,1]<max_stat[,2:ncol(max_stat)])/num_permut

result<-list()
result[[1]]<-p_val
result[[2]]<-max_stat_index
result[[3]]<-critical_value
result[[4]]<-max_stat
path<-paste("iter_",iter,'_type_',type,'.Rdata',sep="")
save(result,file=path)





