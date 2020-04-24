Cluster = function(network, verbose) {
###边聚类
  colnames(network) <- c("node1", "node2")
  x <- network
  # rm(network)
  check.duplicates <- TRUE
  edglim <- 10^4
  if (ncol(x) == 3) {
    wt <- as.numeric(as.character(x[, 3]))
    if (length(which(is.na(wt) == TRUE)) > 0) {
      stop("\nedge weights must be numerical values\n")
    }
    x <- cbind(as.character(x[, 1]), as.character(x[, 2]))
  } else if (ncol(x) == 2) {
    x <- cbind(as.character(x[, 1]), as.character(x[, 2]))
    wt <- NULL
  } else {
    stop("\ninput data must be an edge list with 2 or 3 columns\n")
  }
  if (check.duplicates) { ##去重复边和循环边
    dret <- edge.duplicates(x, verbose = verbose)
    x <- dret$edges
    if (!is.null(wt)) {
      if (length(dret$inds) > 0) {
        wt <- wt[-dret$inds]
      }
    }
    rm(dret)
  }
  el <- x
  len <- nrow(el)
  # nnodes <- length(unique(c(as.character(el[, 1]), as.character(el[,2]))))
  intel <- integer.edgelist(el)
  edges <- intel$edges ##边的数字替代
  node.names <- names(intel$nodes) ##节点名
  numnodes <- length(node.names)
  disk <- FALSE ##未超过边数量，不需要使用磁盘
  dist <- NULL
  directed <- FALSE
  dirweight <- 0.5 
  use.all.edges <- FALSE
  bipartite <- FALSE
  if (is.null(dist)) { ##没有将距离矩阵作为输入
    emptyvec <- rep(1, (len * (len - 1))/2) #全是1
    if (!is.null(wt)) {
      weighted <- TRUE
    }
    else {
      wt <- 0
      weighted <- FALSE
    }
    if (!use.all.edges) { ##仅对有公共端点的两条边计算相似度
      dissvec <- .C("getEdgeSimilarities", as.integer(edges[, 
                                                            1]), as.integer(edges[, 2]), as.integer(len), 
                    rowlen = integer(1), weights = as.double(wt), 
                    as.logical(directed), as.double(dirweight), 
                    as.logical(weighted), as.logical(disk), dissvec = as.double(emptyvec), 
                    as.logical(bipartite), as.logical(verbose))$dissvec
    }
    else { ##对每两条边计算相似度
      dissvec <- .C("getEdgeSimilarities_all", as.integer(edges[, 
                                                                1]), as.integer(edges[, 2]), as.integer(len), 
                    as.integer(numnodes), rowlen = integer(1), 
                    weights = as.double(wt), as.logical(FALSE), 
                    as.double(dirweight), as.logical(weighted), 
                    as.logical(disk), dissvec = as.double(emptyvec), 
                    as.logical(bipartite), as.logical(verbose))$dissvec
    }
    rm(emptyvec)
    distmatrix <- matrix(1, len, len) 
    distmatrix[lower.tri(distmatrix)] <- dissvec  ##将计算的距离转换为下三角矩阵
    rm(dissvec)
    colnames(distmatrix) <- 1:len
    rownames(distmatrix) <- 1:len
    distobj <- as.dist(distmatrix)  ##转换为dist类型的数据
    # rm(distmatrix)
  } else {
    if (class(dist) != "dist") {
      stop("\ndistance matrix must be of class \"dist\" (see ?as.dist)\n")
    }
    else if (attr(dist, which = "Size") != len) {
      stop("\ndistance matrix size must equal the number of edges in the input network\n")
    }
    else if (length(dist) != (len * (len - 1))/2) {
      stop("\ndistance matrix must be the lower triangular matrix of a square matrix\n")
    }
    distobj <- dist
  }
  
  if (verbose) {
    cat("\n")
  }
  rm(distmatrix)
  gc()
  # y <- as.matrix(distobj)+diag(rep(1,len))
  # rm(distobj)
  # gc()
  # y <- 1-y

  return(distobj)
}