
setwd(dirname(rstudioapi::getSourceEditorContext()$path))


source("objTest_fctns.R") 
source("binary_segmentation.R") 
require(sbm)
require(igraph)
require(ecp)


n_nodes<-300
n1=n2=n3=n4=100

m1<-c()
for (i in 1:n1){
  pm <- cbind(c(0.3,0.01,0.01),c(0.01,0.5,0.02),c(0.01,0.02,0.3))
  g1 <- sample_sbm(n_nodes, pref.matrix = pm, block.sizes = c(50,150,100))
  m<-c(as.matrix(as.matrix(g1, "adjacency")))
  m1<-rbind(m1,m)
}
for (i in 1:n2){
  pm <- cbind(c(0.6,0.05,0.02),c(0.05,0.2,0.01),c(0.02,0.01,0.6))
  g2 <- sample_sbm(n_nodes, pref.matrix = pm, block.sizes = c(50,150,100))
  m<-c(as.matrix(as.matrix(g2, "adjacency")))

  m1<-rbind(m1,m)
}

for (i in 1:n3){
  pm <- cbind(c(0.6,0.05,0.02),c(0.05,0.2,0.01),c(0.02,0.01,0.6))
  g3 <- sample_sbm(n_nodes, pref.matrix = pm, block.sizes = c(100,100,100))
  m<-c(as.matrix(as.matrix(g3, "adjacency")))
  m1<-rbind(m1,m)
}
for (i in 1:n4){
  pm <- cbind(c(0.6,0.01),c(0.01,0.3))
  g4 <- sample_sbm(n_nodes, pref.matrix = pm, block.sizes = c(200,100))
  m<-c(as.matrix(as.matrix(g4, "adjacency")))
  m1<-rbind(m1,m)
}

n_nodes=200
m1<-c()
for (i in 1:n1){
  pm <- cbind(c(0.3,0.01,0.01),c(0.01,0.5,0.02),c(0.01,0.02,0.3))
  g1 <- sample_sbm(n_nodes, pref.matrix = pm, block.sizes = c(50,50,100))
  m<-c(as.matrix(as.matrix(g1, "adjacency")))
  m1<-rbind(m1,m)
}

for (i in 1:n2){
  pm <- cbind(c(0.6,0.05,0.02),c(0.05,0.2,0.01),c(0.02,0.01,0.6))
  g2 <- sample_sbm(n_nodes, pref.matrix = pm, block.sizes = c(50,50,100))
  m<-c(as.matrix(as.matrix(g2, "adjacency")))
  
  m1<-rbind(m1,m)
}

for (i in 1:n3){
  pm <- cbind(c(0.6,0.05,0.02),c(0.05,0.2,0.01),c(0.02,0.01,0.6))
  g3 <- sample_sbm(n_nodes, pref.matrix = pm, block.sizes = c(75,75,50))
  m<-c(as.matrix(as.matrix(g3, "adjacency")))
  m1<-rbind(m1,m)
}
for (i in 1:n4){
  pm <- cbind(c(0.6,0.01),c(0.01,0.3))
  g4 <- sample_sbm(n_nodes, pref.matrix = pm, block.sizes = c(150,50))
  m<-c(as.matrix(as.matrix(g4, "adjacency")))
  m1<-rbind(m1,m)
}
save(m1,file="SBM_data_noeds200.Rdata")

distmat<-as.matrix(dist( m1, method = 'manhattan' ))
dim(distmat)
save(distmat,file="SBM_dist_nodes200.Rdata")

y_ecp<-e.divisive(X=m1,sig.lvl=0.01,R=1000,k=NULL,min.size=50)
y_ecp$estimates

#kcp
y_kcp<-kcpa(X=m1, L=6, C=0.05)
y_kcp


iter<-as.numeric(commandArgs(TRUE)[1])
type<-as.numeric(commandArgs(TRUE)[2])

L<-10;minLen<-10
c<-0.1;num_permut<-1000

seed_I<-seeded.intervals(400,decay = sqrt(2))
seed_I<-seed_I[(seed_I[,2]-seed_I[,1])>minLen,]
nrow(seed_I)
#save(seed_I,file="seed_I_SBM.Rdata")
#seed_I<-get(load("seed_I_SBM.Rdata"))


max_stat_index<-c();p_val_interval<-c();critical_value_interval<-c()
type<-1
iter<-200
i<-iter
Data<-m1
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





