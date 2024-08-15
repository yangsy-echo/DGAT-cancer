library(DMwR2)
library(MASS)
library(moments)
library(plyr)
library(parallel)
library(doParallel)
library(foreach)

cl <- makeCluster(45)
registerDoParallel(cl)
cal_randompp <- function(extra_weight,sample_size=25){
  extra_weight$extra_weight <- runif(nrow(extra_weight),min=min(extra_weight$extra_weight),max=max(extra_weight$extra_weight))
  n <- nrow(extra_weight)
  iter <- 0
  dif <- 100
  num0 <- rep(0,n)
  thres <- 1e-4
  ind <- c(1:n)
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
  list(freq)
}
# running example
set.seed(2022)
cancer <- "BLCA"
extra_weight <- read.delim(paste(cancer,"_extra_weight.txt",sep=""),header=T)
#extra_weight$extra_weight <- runif(nrow(extra_weight),min=min(extra_weight$extra_weight),max=max(extra_weight$extra_weight))
gene <- list()
#random_times <- round(1/(0.01/nrow(extra_weight)))
gene <- foreach(x=1:1000000,.combine='c') %dopar% cal_randompp(extra_weight)
gene <- data.frame(gene)
colnames(gene) <- paste("sample",1:1000000,sep = "")
rownames(gene) <- extra_weight$gene
save(gene,file = paste("./random_res/",cancer,"_random.RData",sep = ""))
rm(gene)
stopCluster(cl)

# calculate empirical p values
prop <- read.delim(paste(cancer,"_prop.txt",sep = ""))
load(paste("./random_res/",cancer,"_random.RData",sep = ""))
cat(cancer,"load end.\t")
prop <- prop[match(rownames(gene),prop$gene),]
gene$pp <- prop$pp
gene <- as.matrix(gene)
prop$p.value <- 1
cat("start calculating!\t")
for(j in 1:nrow(gene)){prop$p.value[j] <- sum(gene[j,1:1000000]>=gene[j,1000001])/1000000}
cat("end scoring\n")
#prop$p.value <- apply(gene,1,function(x) sum(x[1:1000000]>=x[1000001])/1000000)
write.table(prop,file = paste(cancer,"_prop_p_value.txt",sep = ""),row.names = F,quote = F,sep = "\t")
rm(gene)
