dunplicate = function(network){
  n <- nrow(network)
  dup <- runif(n,min = 0, max = 0)
  for(i in 1:nrow(network)){
    if(!i%%500){cat("   Checking for loops and duplicate edges ",i*100/n, "%\n")}
    if(network[i,1]==network[i,2]) {
      dup[i] <- 1
      next()
    }
    for(j in i:nrow(network)){
      if(network[i,1]==network[j,2]&&network[i,2]==network[j,1]){
        dup[j] <-1
      }
    }
  }
  return(dup)
}


edge2vec = function(edges,NRvec){
  edge_vec <- data.frame()
  names <- vector()
  for(i in 1:nrow(edges)){
    A <- NRvec[which(rownames(NRvec)==edges[i,1]),]
    B <- NRvec[which(rownames(NRvec)==edges[i,2]),]
    edge_vec<-rbind(edge_vec,cbind(A,B))
    names[i] <- paste(edges[i,1],edges[i,2],sep="__")
  }
  rownames(edge_vec) <- names
  return(edge_vec)
}


density1 = function(network,TF,miRNA,mRNA){
  nnode <- length(TF)+length(miRNA)+length(mRNA)
  TF_ind <- vector()
  for(i in 1:length(TF)){
    TF_ind <- c(TF_ind,union(which(network[,1]==TF[i]),which(network[,2]==TF[i])))
  }
  TF_ind <- sort(TF_ind)
  miRNA_ind <- vector()
  for(i in 1:length(miRNA)){
    miRNA_ind <- c(miRNA_ind,union(which(network[,1]==miRNA[i]),which(network[,2]==miRNA[i])))
  }
  miRNA_ind <- sort(miRNA_ind)
  mRNA_ind <- vector()
  for(i in 1:length(mRNA)){
    mRNA_ind <- c(mRNA_ind,union(which(network[,1]==mRNA[i]),which(network[,2]==mRNA[i])))
  }
  mRNA_ind <- sort(mRNA_ind)
  TF_miRNA_num <- length(intersect(TF_ind,miRNA_ind))
  TF_mRNA_num <- length(intersect(TF_ind,mRNA_ind))
  miRNA_mRNA_num <- length(intersect(miRNA_ind,mRNA_ind))
  mRNA_mRNA_num <- 0
  for(i in mRNA_ind){
    if(length(which(network[i,1]%in%mRNA))!=0&&length(which(network[i,2]%in%mRNA))!=0){
      mRNA_mRNA_num <- mRNA_mRNA_num +1
    }
  }
  num_edge <- TF_miRNA_num + TF_mRNA_num + miRNA_mRNA_num + mRNA_mRNA_num
  
  return(2*num_edge/(nnode*(nnode-1)))
}

Laplace_normalization = function(W){
  DR <- diag(1./sqrt(colSums(W)))
  DR[DR==Inf] <- 0
  DC <- diag(1./sqrt(rowSums(W)))
  DC[DC==Inf] <- 0
  L <- DC %*% W %*% DR
  colnames(L) <- colnames(W)
  rownames(L) <- rownames(W)
  return(L)
}