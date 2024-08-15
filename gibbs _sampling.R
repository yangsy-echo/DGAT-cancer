cal_pp <- function(extra_weight,sample_size=25){
  n <- nrow(extra_weight)
  iter <- 0
  dif <- 100
  num0 <- rep(0,n)
  thres <- 1e-4
  ind <- c(1:n)
  set.seed(123)
  pick <- sample(ind,sample_size,replace = F)
  freq0 <- rep(0,n)
  while(dif>thres){
    iter <- iter+1
    num <- num0
    pick <- sample(ind[-pick],sample_size,replace = F,prob = extra_weight$extra_weight[-pick])
    num[pick] <- num[pick]+1
    if(iter>1) freq0 <- num0/((iter-1)*sample_size)
    freq <- num/(iter*sample_size)
    dif <- norm(freq-freq0,"2")
    num0 <- num
  }
  extra_weight$pp <- freq
  extra_weight
}
# running example
cancer <- "BLCA"
extra_weight <- read.delim(paste(cancer,"_extra_weight.txt",sep=""),header=T,sep="\t")
extra_weight <- cal_pp(extra_weight)
write.table(extra_weight,paste(cancer,"_prop.txt",sep = ""))
