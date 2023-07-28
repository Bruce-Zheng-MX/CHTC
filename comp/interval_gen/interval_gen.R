require("wbs")

require("not")
seeded.intervals <- function(n, decay = sqrt(2), unique.int = F){
  n	<- as.integer(n)
  depth	<- log(n, base = decay)
  depth	<- ceiling(depth)
  
  #n_layer=ceiling(log((min_len-1)/n+1,base=sqrt(0.5)))
  
  
  M	<- sum(2^(1:depth)-1)
  
  boundary_mtx           <- matrix(NA, ncol = 2)
  colnames(boundary_mtx) <- c("st", "end")
  boundary_mtx[1, ]      <- c(1, n)
  
  depth	<- log(n, base = decay)
  depth	<- ceiling(depth)
  
  
  for(i in 2:depth){
    int_length	<- n * (1/decay)^(i-1)
    
    n_int		<- ceiling(round(n/int_length, 14))*2-1		# sometimes very slight numerical inaccuracies
    
    boundary_mtx	<- rbind(boundary_mtx,
                          cbind(floor(seq(1, n-int_length, length.out = (n_int))), 
                                ceiling(seq(int_length, n, length.out = (n_int)))))
  }
  
  if(unique.int){return(unique(boundary_mtx))}
  boundary_mtx
}

random.intervals <- function (n, M, unique.int = F){
  n <- as.integer(n)
  M <- as.integer(M)
  intervals <- matrix(0, nrow = M, ncol = 2)
  intervals[, 1] <- floor(runif(M, min = 1, max = n))		# intended starting point
  intervals[, 2] <- ceiling(runif(M, min = 1, max = n))	# intended end point
  for (i in 1:M){
    tmp <- intervals[i, ] 
    while(tmp[1]==tmp[2]){	                                # in case the two happen to be equal, take new ones
      
      tmp[1] <- floor(runif(1, min = 1, max = n)) 
      tmp[2] <- ceiling(runif(1, min = 1, max = n))
      intervals[i, ] <- tmp
    }
    if(tmp[1] > tmp[2]){intervals[i, ] <- tmp[2:1]}    	# make sure they are ordered
  }
  if(unique.int){return(unique(intervals))}
  intervals
}
seedBS <- function (x, decay = sqrt(2), ...){
  results <- list()
  results$x <- as.numeric(x)
  results$n <- length(results$x)
  results$M <- NA
  results$integrated <- as.logical(FALSE)
  results$rand.intervals <- as.logical(FALSE)
  results$res <- matrix(nrow = 0, ncol = 6)
  if (results$n < 2) 
    stop("x should contain at least two elements")
  if (NA %in% results$x) 
    stop("x vector cannot contain NA's")
  if (var(x) == 0) 
    stop("x is a constant vector, change-point detection is not needed")
  if (decay>2) 
    stop("decay should be <= 2")
  if (decay<=1) 
    stop("decay should be > 1")
  #    if (is.na(results$M)) 
  #        stop("M cannot be NA")
  #    if (length(results$M) > 1) 
  #        stop("M should be a single integer")
  #    if (results$M < 0) 
  #        stop("M should be an integer > 0")
  #    if (results$rand.intervals) 
  #        intervals <- matrix(random.intervals(results$n, results$M), 
  #            ncol = 2)
  #    else {
  intervals <- matrix(seeded.intervals(results$n, decay = decay), 
                      ncol = 2)
  results$M <- nrow(intervals)
  #    }
  #    if (results$integrated) {
  #        results$res <- matrix(.C("wbs_int_rec_wrapper", x = as.double(results$x), 
  #            n = as.integer(results$n), res = double(6 * (results$n - 
  #                1)), intervals = as.integer(intervals), M = as.integer(results$M))$res, 
  #            results$n - 1, 6)
  #    }
  #    else {
  results$res <- matrix(.C("wbs_rec_wrapper", x = as.double(results$x), 
                           n = as.integer(results$n), res = double(6 * (results$n - 
                                                                          1)), intervals = as.integer(intervals), M = as.integer(results$M))$res, 
                        results$n - 1, 6)
  results$res <- matrix(results$res[as.integer(results$res[, 
                                                           1]) > 0, ], ncol = 6)
  #    }
  colnames(results$res) <- c("s", "e", "cpt", "CUSUM", "min.th", 
                             "scale")
  class(results) <- "wbs"
  results$cpt <- changepoints(results)
  return(results)
}

L<-100;minLen<-20
c<-0.1;num_permut<-1000

#seed_I<-seeded.intervals(71-34,decay = sqrt(2))
#seed_I<-seed_I[(seed_I[,2]-seed_I[,1])>minLe

seed_I<-seeded.intervals(170,decay = sqrt(2))
seed_I<-seed_I[(seed_I[,2]-seed_I[,1])>minLen,]
nrow(seed_I)
path<-paste('seed_I_comp_1_170.Rdata')
save(seed_I, file=path)


seed_I<-seeded.intervals(264-171,decay = sqrt(2))
seed_I<-seed_I[(seed_I[,2]-seed_I[,1])>minLen,]
seed_I_2<-seed_I+171
path<-paste('seed_I_172_264.Rdata')
save(seed_I_2, file=path)


#n<-95
#seed_I<-seeded.intervals(n,decay = sqrt(2))
#seed_I<-seed_I[(seed_I[,2]-seed_I[,1])>minLen,]
#path<-paste('seed_I_1_95.Rdata')
#save(seed_I, file=path)

#seed_I<-seeded.intervals(232-95,decay = sqrt(2))
#seed_I<-seed_I[(seed_I[,2]-seed_I[,1])>minLen,]
#seed_I_2<-seed_I+95
#path<-paste('seed_I_96_232.Rdata')
#save(seed_I_2, file=path)

#rand_I<-random.intervals(n,10000)
#rand_I<-rand_I[(rand_I[,2]-rand_I[,1])>minLen,][1:L,]
#rand_I<-rbind(rand_I,c(1,232))

#path<-paste('seed_I_1_95.Rdata')
#save(seed_I, file=path)
#path<-paste('rand_I.Rdata')
#save(rand_I, file=path)
