train = function(tf_miRNA,tf_mRNA,miRNA_mRNA,mRNA_mRNA,miRNA,TF,mRNA,dim,k_max){
  # setwd("F:\\实验室\\模块\\论文二\\GUI\\shiny-download")
  # library(plyr)
  # library(e1071)
  # library(pROC)
  # source("function.R")
  # source("train.R")
  # 
  # tf_miRNA <-read.csv("tf_miRNA_small.txt", header = F, sep = "\t", stringsAsFactors = F)
  # tf_mRNA <- read.csv("tf_mRNA_small.txt", header = F, sep = "\t", stringsAsFactors = F)
  # mRNA_mRNA <- read.csv("mRNA_mRNA_small.txt", header = F, sep = "\t", stringsAsFactors = F)
  # miRNA_mRNA <- read.csv("miRNA_mRNA_small.txt", header = F, sep = "\t", stringsAsFactors = F)
  # 
  # miRNA <- read.csv("BRCA_miRNA.csv",header = TRUE, row.names = 1, sep = ",",stringsAsFactors = F)
  # TF <-read.csv("BRCA_tf.csv",header = TRUE, row.names = 1, sep = ",",stringsAsFactors = F)
  # mRNA <- read.csv("BRCA_mRNA.csv",header = TRUE, row.names = 1, sep = ",",stringsAsFactors = F)
  
  cat("   Train\n")
  
  k<-5
  nr <- nrow(miRNA_mRNA)
  cvlist <- list()
  set.seed(10)
  n <- rep(1:k,ceiling(nr/k))[1:nr]
  temp <- sample(n,nr)
  x <- 1:k
  dataseq <- 1:nr
  cvlist <- lapply(x, function(x) dataseq[temp==x])

  vec_filename <- list()
  train_id <-  list()
  test_id <- list()
  filename<-list()
  dense_filename <-list()
  AUC <- vector()

  for(p in 1:length(k_max)){
    for(q in 1:length(dim)){
      for(i in 1:k){
        test_id[[i]] <- cvlist[[i]]
        train_id[[i]] <- setdiff(c(1:nr),test_id[[i]])
        filename[[i]] <- paste("network_train",i,".txt",sep = "")
        dense_filename[[i]] <- paste("network_train",i,"_dense.txt",sep = "")
        vec_filename[[i]] <- paste("network_train",i,"_vec.txt",sep = "")
      }

      for(i in 1:k){
        test_link <- miRNA_mRNA[test_id[[i]],]
        train_link <- miRNA_mRNA[train_id[[i]],]
        network <- rbind(train_link,tf_miRNA,tf_mRNA,mRNA_mRNA)
        n<-network[,c(2,1)]
        # colnames(n)<-colnames(network)
        network_end <- rbind(network,n)
        network_end <- cbind(network_end,runif(nrow(network_end),min = 1, max = 1))
        write.table(network_end, file = filename[[i]], sep = "\t", row.names = FALSE,
                    col.names = FALSE, quote = FALSE)
        command_1 <- paste("reconstruct.exe -train ",filename[[i]]," -output ",
                           dense_filename[[i]]," -depth 2 -threshold ",k_max[p],sep = "")
        system(command = command_1, intern = FALSE, wait = TRUE, invisible = TRUE)
        command_2 <- paste("line.exe -train ",dense_filename[[i]]," -output ",
                           vec_filename[[i]]," -binary 0 -size ",dim[q],"-order 2
                         -negative 5 -samples 100 -threads 10",sep = "")
        system(command = command_2, intern = FALSE, wait = TRUE, invisible = TRUE)
      }

      negative_ind <- list()

      for(i in 1:k){
        vec_filename[[i]] <- paste("network_train",i,"_vec.txt",sep = "")
        NRvec <- read.csv(vec_filename[[i]], header = F, sep = " ", stringsAsFactors = F)
        NRvec <- NRvec[-1,]
        rownames(NRvec) <- NRvec[,1]
        NRvec <- NRvec[,-c(1,ncol(NRvec))]
        net_miname <- vector()
        net_mname <- vector()
        for(j in 1:nrow(NRvec)){
          if(rownames(NRvec)[j]%in%rownames(miRNA)) net_miname <- rbind(net_miname,rownames(NRvec)[j])
          if(rownames(NRvec)[j]%in%rownames(mRNA)) net_mname <- rbind(net_mname,rownames(NRvec)[j])
        }
        net_miname <- unique(as.vector(net_miname))
        net_mname <- unique(as.vector(net_mname))
        noedges <- data.frame()
        for(j in 1:length(net_miname)){
          for(g in 1:length(net_mname)){
            if(length(intersect(which(tf_miRNA[,2]==net_miname[j]),which(tf_miRNA[,1]==net_mname[g])))==0){
              noedges <- rbind(noedges,cbind(net_miname[j],net_mname[g]))
            }
          }
        }
        noed_vec <- edge2vec(noedges,NRvec)
        num_negative <- nrow(miRNA_mRNA)-length(cvlist[[i]])
        negative_ind[[i]] <- sample(1:nrow(noedges), num_negative)
        negative_vec <- noed_vec[negative_ind[[i]],]
        negative_vec <- cbind(negative_vec,runif(nrow(negative_vec),min = 0,max = 0)) ##训练负样本
        colnames(negative_vec)[ncol(negative_vec)]<-"label"
        positive_edges <- miRNA_mRNA[train_id[[i]],]
        positive_vec <- edge2vec(positive_edges,NRvec)
        positive_vec1 <- cbind(positive_vec,runif(nrow(positive_vec),min = 1,max = 1)) ##训练正样本
        colnames(positive_vec1)[ncol(positive_vec1)]<-"label"
        colnames(negative_vec) <- colnames(positive_vec1)
        train <- rbind(negative_vec,positive_vec1)  ##训练样本
        test_positive <- miRNA_mRNA[test_id[[i]],]
        positive_test <- data.frame()
        for(j in 1:nrow(test_positive)){
          if(length(which(rownames(NRvec)==test_positive[j,1]))!=0&&
             length(which(rownames(NRvec)==test_positive[j,2]))!=0){
            positive_test <- rbind(positive_test,test_positive[j,])
          }
        }
        positive_test_vec <- edge2vec(positive_test,NRvec)
        colnames(positive_test_vec) <- colnames(positive_vec)
        positive_test_vec <- rbind(positive_test_vec,positive_vec)  ##测试正样本
        negative_test_vec <- noed_vec ##测试负样本
        test_vec <- rbind(positive_test_vec,negative_test_vec)
        colnames(test_vec) <- c(1:ncol(test_vec))
        rownames(test_vec) <- c(1:nrow(test_vec))
        test_labels <- as.vector(c(rep(1,times = nrow(positive_test_vec)),rep(0,times = nrow(negative_test_vec))))
        colnames(train)[1:(ncol(train)-1)] <- c(1:(ncol(train)-1))
        train_labels <- train$label
        train <- train[,-ncol(train)]
        rownames(train) <- c(1:nrow(train))
        ##预测
        model<- svm(train,train_labels,kernel = 'radial', gamma = if (is.vector(train)) 1 else 1 / ncol(train) )
        pre_svm <- predict(model,test_vec)

        ROC <- roc(test_labels,pre_svm)
        AUC[i] <- ROC$auc
      }
      # plot(ROC,print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),
      #      grid.col=c("green", "red"), max.auc.polygon=TRUE,
      #      auc.polygon.col="skyblue", print.thres=TRUE)
      avg_auc <- sum(AUC)/k
      # cat("dim =",dim, "k_max =",k_max, "avg_auc =", avg_auc,"\n")
    }

  }
  cat("   Train over\n")
  return(ROC)
}