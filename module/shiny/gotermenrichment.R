goanalysis = function(targetgene,method){
  # setwd( "F:\\实验室\\模块\\论文二\\GUI\\shiny-download")
  # 
  ## 读入数据
  if(is.null(targetgene)) {
    if(method == "LinkModule"){
      targetgene <- read.table("targetgene.txt",header = F, sep = "\t")
    } else if(method == "LinkPartition"){
      targetgene <- read.table("targetgene_LinkPartition.txt",header = F, sep = "\t")
    } else if(method == "DBSCAN"){
      targetgene <- read.table("targetgene_DBSCAN.txt",header = F, sep = "\t")
    }
  }
  # if(is.null(targetgene)){
  #   targetgene <- read.table("targetgene.txt",header = F, sep = "\t")
  # }
  
  # ## 加载包
  # library("org.Hs.eg.db")
  # library("GSEABase")
  # library("GOstats")
  
  ## 初始化表头
  goAnn <- get("org.Hs.egGO")
  universe <- Lkeys(goAnn)
  GOcolnames <- c("GOBPID","Pvalue","OddsRatio","ExpCount","Count","Size","Term","Symbols","GeneRatio")
  gorownames <- c("module",GOcolnames)


  ## 对target gene进行富集分析
  # geneGO <- matrix(,,length(gorownames))
  # colnames(geneGO) <- gorownames
  # geneGO <- geneGO[-1,]
  # 
  # 
  # for(i in 1:length(targetgene[,1])){
  #   #转成字符格式，忽略因子
  #   temp <- as.character(targetgene[i,1])
  #   #根据逗号分割出全部基因，变成list，整体作为一个元素放在temp的第一个元素中
  #   temp <- strsplit(temp,split = " ")
  # 
  #   gene <- temp[[1]]
  # 
  #   if(gene!=""){
  #     entrezIDs <- mget(gene, org.Hs.egSYMBOL2EG, ifnotfound=NA)
  #     entrezIDs <- as.character(entrezIDs)
  # 
  #     ##GO BP enrichment analysis
  #     params <- new("GOHyperGParams",
  #                   geneIds=entrezIDs,
  #                   universeGeneIds=universe,
  #                   annotation="org.Hs.eg.db",
  #                   ontology="BP",
  #                   pvalueCutoff=0.05,
  #                   conditional=FALSE,
  #                   testDirection="over")
  # 
  #     over1 <- hyperGTest(params)
  #     glist <- geneIdsByCategory(over1)
  #     glist <- sapply(glist, function(.ids) {
  #       .sym <- mget(.ids, envir=org.Hs.egSYMBOL, ifnotfound=NA)
  #       .sym[is.na(.sym)] <- .ids[is.na(.sym)]
  #       paste(.sym, collapse=";")
  #     })
  # 
  #     go <-summary(over1)
  #     go$Symbols <- glist[as.character(go$GOBPID)] #富集基因
  # 
  #     go <- go[go$Count>1,] #count小于2的去掉
  # 
  #     go$GeneRatio <- go$Count/go$Size
  # 
  #     go <- as.matrix(go)
  # 
  #     if(class(go)=="character")
  #       go <- t(as.matrix(go))
  # 
  #   }else{
  #     go <- matrix(0,nrow=1,ncol=9)
  #     colnames(go)=GOcolnames
  #   }
  # 
  #   if(!nrow(go)){
  #     go <-matrix(0,nrow=1,ncol=9)
  #     colnames(go)=GOcolnames
  #   }
  # 
  #   tGO<-cBind(i,go)
  #   cat("i: ", i, "\n")
  #   colnames(tGO) <- gorownames
  #   geneGO <- rBind(geneGO,tGO)
  # }
  # write.csv(geneGO, file="GO.csv")
  # return(geneGO)
  # #输出文件
  # write.csv(geneGO, file="BRCALinkModuleGOterms.csv")
  
  
  
  ######################################################
  ### KEGG  ##################################
  ##################################
  ######################
  ##############
  ########
  KEGGAnn <- get("org.Hs.egPATH")
  universe <- Lkeys(KEGGAnn)
  KEGGcolnames <- c("KEGGID","Pvalue","OddsRatio","ExpCount","Count","Size","pathway","Symbols","GeneRatio")
  KEGGrownames <- c("module",KEGGcolnames)


  ## 对target gene进行富集分析
  geneKEGG <- matrix(,,length(KEGGrownames))
  colnames(geneKEGG) <- KEGGrownames
  geneKEGG <- geneKEGG[-1,]


  for(i in 1:length(targetgene[,1])){
    #转成字符格式，忽略因子
    temp <- as.character(targetgene[i,1])
    #根据逗号分割出全部基因，变成list，整体作为一个元素放在temp的第一个元素中
    temp <- strsplit(temp,split = " ")

    gene <- temp[[1]]

    if(gene!=""){
      entrezIDs <- mget(gene, org.Hs.egSYMBOL2EG, ifnotfound=NA)
      entrezIDs <- as.character(entrezIDs)

      ##KEGG enrichment analysis
      params <- new("KEGGHyperGParams",
                    geneIds=entrezIDs,
                    universeGeneIds=universe,
                    annotation="org.Hs.eg.db",
                    categoryName="KEGG",
                    pvalueCutoff=0.05,
                    testDirection="over")
      over <- hyperGTest(params)
      kegg <- summary(over)
      library(Category)
      glist <- geneIdsByCategory(over)
      glist <- sapply(glist, function(.ids) {
        .sym <- mget(.ids, envir=org.Hs.egSYMBOL, ifnotfound=NA)
        .sym[is.na(.sym)] <- .ids[is.na(.sym)]
        paste(.sym, collapse=";")
      })
      kegg$Symbols <- glist[as.character(kegg$KEGGID)]

      kegg <- kegg[kegg$Count>1,] #count小于2的去掉

      kegg$GeneRatio <- kegg$Count/kegg$Size

      kegg <- as.matrix(kegg)

      if(class(kegg)=="character")
        kegg <- t(as.matrix(kegg))

    }else{
      kegg <- matrix(0,nrow=1,ncol=9)
      colnames(kegg) = KEGGcolnames
    }

    if(!nrow(kegg)){
      kegg <- matrix(0,nrow=1,ncol=9)
      colnames(kegg) = KEGGcolnames
    }

    tKEGG <- cBind(i, kegg)
    cat("i: ", i, "\n")
    colnames(tKEGG) <- KEGGrownames
    geneKEGG <- rBind(geneKEGG, tKEGG)
  }
  write.csv(geneKEGG, file="KEGG.csv")
  return(geneKEGG)
  # write.csv(geneKEGG, file="BRCALinkModuleKEGG.csv")
}


