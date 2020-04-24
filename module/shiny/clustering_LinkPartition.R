LinkPartition = function(tf_miRNA,tf_mRNA,miRNA_mRNA,mRNA_mRNA,miRNA,TF,mRNA,network_similarity){
  # setwd("F:\\实验室\\模块\\论文二\\GUI\\shiny-download")
  # ###LinkPartition
  # library(plyr)
  # library(linkcomm)
  # source("function.R")
  # source("LinkCluster.R")
  cat("   LinkPartition\n")
  # if(is.null(tf_miRNA)||is.na(tf_mRNA)||is.na(miRNA_mRNA)||is.na(mRNA_mRNA)||is.na(miRNA)||is.na(TF)||is.na(mRNA)||is.na(network_similarity)){
  #   if(file.exists("cluster.result_LinkPartition.txt")){
  #     file.remove("cluster.result_LinkPartition.txt")
  #   }
  #   sink("cluster.result_LinkPartition.txt",append = TRUE, split = FALSE)
  #   cat("NULL")
  #   sink()
  #   if(file.exists("targetgene_LinkPartition.txt")){
  #     file.remove("targetgene_LinkPartition.txt")
  #   }
  #   sink("targetgene_LinkPartition.txt",append = TRUE, split = FALSE)
  #   cat("NULL")
  #   sink()
  #   result=NULL
  #   cat("   LinkPartition don't have input\n")
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
  
  network <- rbind(tf_miRNA,miRNA_mRNA,tf_mRNA,mRNA_mRNA)
  verbose <- TRUE
  cluster <- getLinkCommunities(network,plot = FALSE)
  clus <- cluster$clusters
  # distobj <- Cluster(network, verbose)
  # hcedges <- hclust(distobj, method = "average")
  # cluster <- hcedges$labels
  
  if (verbose) {
    cat("   Processing\n")
  }
  # lclus <- unique(cluster)
  # clus <- list()
  # for (i in 1:length(lclus)) {
  #   clus[[i]] <- which(cluster == lclus[i])
  # }
  
  clusters <- list()
  for (i in 1:length(clus)) {
    names <- vector()
    for(j in 1:length(clus[[i]])){
      a <- c(network[clus[[i]][j],1],network[clus[[i]][j],2])
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
  # cat("avg_density = ", mean(cludensity))
  
  if(file.exists("cluster.result_LinkPartition.txt")){
    file.remove("cluster.result_LinkPartition.txt")
  }
  names(last_cluster) <- paste("M",1:length(last_cluster),sep="")
  sink("cluster.result_LinkPartition.txt",append = TRUE, split = FALSE)
  for(i in 1:length(last_cluster)){
    cat(names(last_cluster)[i],"  density = ",cludensity[i],"\n")
    cat("TF    ",length(TF_clu[[i]]), "    ",TF_clu[[i]],"\n")
    cat("miRNA    ",length(miRNA_clu[[i]]), "    ",miRNA_clu[[i]],"\n")
    cat("mRNA    ",length(mRNA_clu[[i]]), "    ",mRNA_clu[[i]],"\n\n")
  }
  cat("avg_density = ", mean(cludensity))
  sink()
  
  if(file.exists("targetgene_LinkPartition.txt")){
    file.remove("targetgene_LinkPartition.txt")
  }
  sink("targetgene_LinkPartition.txt",append = TRUE, split = FALSE)
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
  cat("   LinkPartition over\n")
  return(result)
}