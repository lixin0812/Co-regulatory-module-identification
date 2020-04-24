dbscancluster = function(tf_miRNA,tf_mRNA,miRNA_mRNA,mRNA_mRNA,miRNA,TF,mRNA,network_similarity){
  # setwd("F:\\实验室\\模块\\论文二\\GUI\\shiny-download")
  ###DBSCAN+降维+加权求和+非边特征聚类
  # library(plyr)
  # library(dbscan)
  # source("function.R")
  cat("   DBSCAN\n")
  # if(is.null(tf_miRNA)||is.na(tf_mRNA)||is.na(miRNA_mRNA)||is.na(mRNA_mRNA)||is.na(miRNA)||is.na(TF)||is.na(mRNA)||is.na(network_similarity)){
  #   if(file.exists("cluster.result_DBSCAN.txt")){
  #     file.remove("cluster.result_DBSCAN.txt")
  #   }
  #   sink("cluster.result_DBSCAN.txt",append = TRUE, split = FALSE)
  #   cat("NULL")
  #   sink()
  #   if(file.exists("targetgene_DBSCAN.txt")){
  #     file.remove("targetgene_DBSCAN.txt")
  #   }
  #   sink("targetgene_DBSCAN.txt",append = TRUE, split = FALSE)
  #   cat("NULL")
  #   sink()
  #   result=NULL
  #   cat("   DBSCAN don't have input\n")
  #   return(result)
  # }
  # tf_miRNA_BRCA <-read.csv("tf_miRNA_BRCA.txt", header = F, sep = "\t", stringsAsFactors = F)
  # tf_mRNA_BRCA <- read.csv("tf_mRNA_BRCA.txt", header = F, sep = "\t", stringsAsFactors = F)
  # mRNA_mRNA_BRCA <- read.csv("mRNA_mRNA_BRCA.txt", header = F, sep = "\t", stringsAsFactors = F)
  # miRNA_mRNA_BRCA <- read.csv("miRNA_mRNA_BRCA.txt", header = F, sep = "\t", stringsAsFactors = F)
  # 
  # miRNA <- read.csv("BRCA_miRNA.csv",header = TRUE, row.names = 1, sep = ",",stringsAsFactors = F)
  # TF <-read.csv("BRCA_tf.csv",header = TRUE, row.names = 1, sep = ",",stringsAsFactors = F)
  # mRNA <- read.csv("BRCA_mRNA.csv",header = TRUE, row.names = 1, sep = ",",stringsAsFactors = F)
  # 
  # network_similarity <- read.csv(file = "BRCA_gene_MIC_lfdr.csv", header = T, sep = ",", row.names = 1)
  colnames(network_similarity) <- rownames(network_similarity)
  
  network <- rbind(tf_miRNA,miRNA_mRNA,tf_mRNA,mRNA_mRNA)
  network<-cbind(network, runif(nrow(network),min = 1, max = 1))
  network2 <- network[,c(2,1,3)]
  colnames(network2) <- colnames(network)
  network_all <- rbind(network,network2)
  
  write.table(network_all, file = "network_DBSCAN.txt", sep = "\t", row.names = FALSE,
              col.names = FALSE, quote = FALSE)
  vec_filename <- list()
  for(i in 1:5){
    vec_filename[[i]] <- paste("NRvec_DBSCAN_",i,".txt",sep = "")
  }
  command_1 <- "reconstruct.exe -train network_DBSCAN.txt -output network_dense_DBSCAN.txt -depth 2 -threshold 50"
  system(command = command_1, intern = FALSE, wait = TRUE, invisible = TRUE)
  for(i in 1:5){
    command_2 <- paste("line.exe -train network_dense_DBSCAN.txt -output ",
                       vec_filename[[i]],
                       " -binary 0 -size 80 -order 2 -negative 5 -samples 100 -threads 5",
                       sep = "")
    system(command = command_2, intern = FALSE, wait = TRUE, invisible = TRUE)
  }
  for(i in 1:5){
    NRvec_BRCA1 <- read.csv(vec_filename[[i]], header = F, sep = " ", stringsAsFactors = F)
    NRvec_BRCA1 <- NRvec_BRCA1[-1,]
    rownames(NRvec_BRCA1) <- NRvec_BRCA1[,1]
    NRvec_BRCA1 <- NRvec_BRCA1[,-c(1,ncol(NRvec_BRCA1))]
    if(i==1){
      NRvec_BRCA2 <- NRvec_BRCA1
    } else {
      NRvec_BRCA2 <- NRvec_BRCA2 + NRvec_BRCA1
    }
  }
  NRvec_BRCA <- NRvec_BRCA2/5
  write.table(NRvec_BRCA, file = "NRvec_DBSCAN_avg.txt", sep = " ", row.names = TRUE,
              col.names = FALSE, quote = FALSE)
  
  
  NRvec_BRCA <- read.csv("NRvec_DBSCAN_avg.txt", header = F,row.names = 1, sep = " ", stringsAsFactors = F)
  gene_names<-unique(c(rownames(NRvec_BRCA),network[,1],network[,2]))
  a<-which(rownames(network_similarity)%in%gene_names)
  similarity <- network_similarity[a,]
  similarity <- similarity[,a]
  NRvec_BRCA <- NRvec_BRCA[rownames(similarity),]
  NRvec_BRCA1 <- as.data.frame(scale(NRvec_BRCA))
  
  
  ##降维
  similarity1 <- as.data.frame(cmdscale(similarity,k=80))
  similarity1 <- as.data.frame(scale(similarity1))
  verbose <- TRUE
  if (verbose) {
    cat("   edge2vec\n")
  }
  ws <- 0.5
  feature <- ws*similarity1 + (1-ws)*NRvec_BRCA1
  # feature <- cbind(similarity1,NRvec_BRCA1)
  # network_vec <- edge2vec(network_all,feature)
  network_vec <- feature
  
  verbose <- TRUE
  if (verbose) {
    cat("   Clustering\n")
  }
  # if(file.exists("density_DBSCAN_withoutlink.txt")){
  #   file.remove("density_DBSCAN_withoutlink.txt")
  # }
  knns <- 5
  kNNdistplot(network_vec, k = knns)
  dist <- sort(kNNdist(network_vec, knns))
  # dist <- sort(kNNdist(network_vec, 5),decreasing = TRUE)
  dist_mean <- mean(dist)
  # plot(sort(kNNdist,decreasing = TRUE), type = "l", ylab = paste(5, "-NN distance",
  #                                              sep = ""), xlab = "Points (sample) sorted by distance")
  # abline(h = dist_mean, lty = 3)
  # g<-0
  # for(i in 1:(length(dist)-1)){
  #   if(as.numeric(dist[i+1])==as.numeric(dist[i])){
  #     next()
  #   }
  #   temp <- dist[i+1]-dist[i]
  #   if(temp>=g){
  #     g<-temp
  #     min<-i
  #   }
  #   # } else{
  #   #   eps <- kNNdist[i]
  #   #   break
  #   # }
  # }
  # if(min>dist_mean){
  #   epss <- dist_mean
  # } else epss <- dist[min]
  epss <- dist_mean
  # for(minpt in 5){
  #   for(eps in seq(1,15,0.1)){
  for(minpt in 5){
    for(eps in 4){
      cluster.out <- dbscan(network_vec, eps = eps, minPts = minpt)
      cluster <- cluster.out$cluster
      
      
      # if (verbose) {
      #   cat("   Processing\n")
      # }
      lclus <- unique(cluster)
      # lclus <- lclus[which(lclus>0)]
      # if(length(lclus)==0) next
      clus <- list()
      for (i in 1:length(lclus)) {
        clus[[i]] <- which(cluster == lclus[i])
      }
      
      clusters <- list()
      for (i in 1:length(clus)) {
        names <- vector()
        for(j in 1:length(clus[[i]])){
          a <- c(network_all[clus[[i]][j],1],network_all[clus[[i]][j],2])
          names <- c(names,a)
        }
        clusters[[i]] <- unique(names)
      }
      last_cluster<-list()
      num <- 1
      for(i in 1:length(clusters)){
        mi <- FALSE
        tf <- FALSE
        m <- FALSE
        for(j in 1:length(clusters[[i]])){
          if(tf&&mi&&m){
            break
          } else if(!tf&&(clusters[[i]][j] %in% rownames(TF))){
            tf <- TRUE
          } else if(!mi&&(clusters[[i]][j] %in% rownames(miRNA))){
            mi <- TRUE
          } else if(!m&&((clusters[[i]][j] %in% rownames(mRNA)))){
            m <-TRUE
          }
        }
        if(tf&&mi&&m){
          last_cluster[[num]] <- clusters[[i]]
          num <- num+1
        }
      }
      if(length(last_cluster)==0) next
      for(p in 1:length(last_cluster)){
        q <- p+1
        while(q <= length(last_cluster)){
          un <- length(union(last_cluster[[p]],last_cluster[[q]]))
          int <- length(intersect(last_cluster[[p]],last_cluster[[q]]))
          if(un ==int){
            last_cluster[[q]] <- NULL
          } else {
            q <- q+1
          }
        }
      }
      
      TF_clu <- list()
      miRNA_clu<-list()
      mRNA_clu <- list()
      for(i in 1:length(last_cluster)){
        TF_clu[[i]] <- last_cluster[[i]][which(last_cluster[[i]] %in% rownames(TF))]
        miRNA_clu[[i]] <- last_cluster[[i]][which(last_cluster[[i]] %in% rownames(miRNA))]
        mRNA_clu[[i]] <- last_cluster[[i]][which(last_cluster[[i]] %in% rownames(mRNA))]
      }
      cludensity <- vector()
      for(i in 1: length(last_cluster)){
        cludensity[i] <- density1(network,TF_clu[[i]],miRNA_clu[[i]],mRNA_clu[[i]])
      }
      cat("minpt = ",minpt,"eps = ",eps,"avg_density = ", mean(cludensity),"\n")
      # sink("density_DBSCAN_withoutlink.txt",append = TRUE, split = FALSE)
      # cat("minpt = ",minpt,"eps = ",eps,"avg_density = ", mean(cludensity),"\n")
      # sink()
    }
  }
  
  if(file.exists("cluster.result_DBSCAN.txt")){
    file.remove("cluster.result_DBSCAN.txt")
  }
  names(last_cluster) <- paste("M",1:length(last_cluster),sep="")
  sink("cluster.result_DBSCAN.txt",append = TRUE, split = FALSE)
  for(i in 1:length(last_cluster)){
    cat(names(last_cluster)[i],"  density = ",cludensity[i],"\n")
    cat("TF    ",length(TF_clu[[i]]), "    ",TF_clu[[i]],"\n")
    cat("miRNA    ",length(miRNA_clu[[i]]), "    ",miRNA_clu[[i]],"\n")
    cat("mRNA    ",length(mRNA_clu[[i]]), "    ",mRNA_clu[[i]],"\n")
  }
  cat("avg_density = ", mean(cludensity))
  sink()
  
  if(file.exists("targetgene_DBSCAN.txt")){
    file.remove("targetgene_DBSCAN.txt")
  }
  sink("targetgene_DBSCAN.txt",append = TRUE, split = FALSE)
  for(i in 1:length(last_cluster)){
    if(cludensity[i]<0.015) next
    cat(mRNA_clu[[i]],"\n")
  }
  sink()
  result <- data.frame()
  ro <- 1
  for (i in 1:length(last_cluster)) {
    temp <- paste(names(last_cluster)[i],"  density = ",cludensity[i],sep = "")
    result[ro,1] <- temp
    ro <- ro+1
    TF_head <- paste("TF    ",length(TF_clu[[i]])," ",sep = "")
    TF_name <- vector()
    for (j in 1:length(TF_clu[[i]])) {
      TF_name <- paste(TF_name," ",TF_clu[[i]][j],sep = "")
    }
    TF_name <- paste(TF_head,TF_name,sep = "    ")
    result[ro,1] <- TF_name
    ro <- ro+1
    miRNA_head <- paste("miRNA    ",length(miRNA_clu[[i]])," ",sep = "")
    miRNA_name <- vector()
    for (j in 1:length(miRNA_clu[[i]])) {
      miRNA_name <- paste(miRNA_name," ",miRNA_clu[[i]][j],sep = "")
    }
    miRNA_name <- paste(miRNA_head,miRNA_name,sep = "    ")
    result[ro,1] <- miRNA_name
    ro <- ro+1
    mRNA_head <- paste("mRNA    ",length(mRNA_clu[[i]])," ",sep = "")
    mRNA_name <- vector()
    for (j in 1:length(mRNA_clu[[i]])) {
      mRNA_name <- paste(mRNA_name," ",mRNA_clu[[i]][j],sep = "")
    }
    mRNA_name <- paste(mRNA_head,mRNA_name,sep = "    ")
    result[ro,1] <- mRNA_name
    ro <- ro+1
  }
  result[ro,1] <- paste("avg_density = ", mean(cludensity))
  cat("   DBSCAN over\n")
  return(result)
}