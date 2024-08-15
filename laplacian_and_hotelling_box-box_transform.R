library(DMwR2)
library(MASS)
library(moments)

score_func <- function(score_df,L_m=20,m){
  # L_m was the number of selected most important features in laplacian selection
  # m was the number of chosen principal components
  cal_num <- unlist(apply(score_df,1,function(x){sum(is.na(x))}))
  score_df2 <- score_df[cal_num <= ncol(score_df)/2,]
  score_df2$JSD <- 1-sqrt(score_df2$JSD)
  score_df2$C_score <- 1-sqrt(score_df2$C_score)
  score_mat <- scale(score_df2,center = T,scale = T)
  # knn imputation
  knn.full <- knnImputation(score_mat,k = 5,scale = F,meth = "weighAvg")
  # k nearest neighbors
  dist.mat <- dist(knn.full,upper = F,p=2)
  dist.mat <- as.matrix(dist.mat)
  #k <- 15
  k <- round(0.01*nrow(score_df2))
  knn.list <- list()
  gene_name <- data.frame(gene=rownames(dist.mat),ind=1:nrow(dist.mat))
  for(i in 1:nrow(dist.mat)){
    knn.list[[i]] <- gene_name$ind[order(dist.mat[,i],decreasing = F)[2:(1+k)]]
  }
  names(knn.list) <- gene_name$gene
  # construct weight network
  S <- matrix(0,nrow = nrow(score_df2),ncol = nrow(score_df2))
  tt <- 1
  for(i in 1:nrow(score_df2)){
    S[i,knn.list[[i]]] <-S[knn.list[[i]],i] <- exp(-dist.mat[i,knn.list[[i]]]^2/tt);
  }
  # compute score
  Fea <- knn.full
  Dd <- diag(rowSums(S))
  L <- Dd - S
  fr_hat <- function(x){
    mm <- (t(x)%*%rowSums(Dd))/sum(Dd)
    x-rep(mm,length(x))
  }
  fr <- apply(Fea,2,FUN = fr_hat)
  Lr_func <- function(x){
    (t(x)%*%L%*%x)/(t(x)%*%Dd%*%x)
  }
  Lr <- apply(fr,2,FUN = Lr_func)
  selected2 <- colnames(score_df)[order(Lr,decreasing = F)[1:L_m]]
  cal_num <- unlist(apply(score_df[,selected2],1,function(x){sum(is.na(x))}))
  score_df2 <- score_df[cal_num <= length(selected2)/2,selected2]
  
  if("JSD"%in%colnames(score_df2)){score_df2$JSD <- 1-sqrt(score_df2$JSD)}
  if("C_score"%in%colnames(score_df2)){score_df2$C_score <- 1-sqrt(score_df2$C_score)}
  score_mat <- scale(score_df2,center = T,scale = T)
  # knn imputation
  knn.full <- knnImputation(score_mat,k = 5,scale = F,meth = "weighAvg")
  data <- knn.full
  C <- cov(data)
  M <- colMeans(data)
  a.e <- eigen(C,symmetric = T)
  V <- a.e$vectors
  #result <- matrix(nrow = nrow(V),ncol = 20)
  #cat(m,"\n")
  U <- V[1:m,]
  data_h <- U%*%(apply(data,1,function(x) x-M))
  data_h <- t(data_h)
  # box-cox transform
  for(i in 1:ncol(data_h)){
    lam <- -floor(min(c(data_h[,i],0)))
    y <- data_h[,i] + lam
    b <- boxcox(y~1)
    lambda <- b$x[which.max(b$y)]
    data_h[,i] <- y^lambda
  }
  data_h <- scale(data_h)
  data_p <- apply(data_h,2,function(x) unlist(sapply(x,pnorm)))
  
  extra_weight<-(-log(apply(data_p,1,function(x) { pchisq(-2*sum(log(x)),df=2*length(x),lower.tail=F)})))
  extra_weight <- data.frame(gene=rownames(score_df2),extra_weight=extra_weight)
  extra_weight <- extra_weight[match(unique(extra_weight$gene),extra_weight$gene),]
  extra_weight
}
# running example
extra_weight <- score_func(cancer.list[[1]],m=10)
cancer <- "BLCA"
write.table(extra_weight,file=paste(cancer,"_extra_weight.txt",sep=""),sep="\t",row.names=F)
