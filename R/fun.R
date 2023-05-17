#' @import Seurat
#' @import magrittr
#' @importFrom pbmcapply pbmclapply
#' @importFrom Matrix rowMeans
#' @importFrom stats aggregate as.dist kmeans setNames
#' @importFrom utils head
#' @importFrom rlang .data
#' @include AllGenerics.R utils.R
#' @noRd
NULL


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#'
#' @param assay character.\cr Assay of Seurat object used to calculate sub-cluster
#'   expression matrix. Default is "RNA".
#' @param reduction character.\cr Reduction of Seurat object used to calculate
#'   reduciont position of sub-cluster. Default is "umap".
#' @param cell.range integer.\cr Range of change of cells number in a sub-cluster
#'   when runing Kmean. Default is 20.
#' @param cell.size integer.\cr Number of cells in a sub-cluster when runing Kmean.
#'   Default is 50.
#' @param set.seed integer.\cr Seed number for Kmean to compute sub-cluster of cell cluster.
#'   Default is 1.
#' @param n.core integer.\cr Number of core for parallel job of running MinMax Kmean cluster.
#'   Default is 1.
#'
#' @rdname CreateSubclusterObject
#' @export
#' @method CreateSubclusterObject Seurat
#'
CreateSubclusterObject.Seurat <- function(object, assay = "RNA", reduction = "umap", cell.range = 20, cell.size = 50, set.seed = 1, n.core = 1, ...){
  ## runSeuratWorkflow #####
  object <- runSeuratWorkflow(object = object, assay = assay)
  if(reduction == "umap"){
    if(!"umap" %in% names(object@reductions))
      object <- RunUMAP(object, reduction = "pca", dims = 1:30, verbose = F)
  }else{
    # if(!reduction %in% object@reductions) stop(reduction, " not in SeuratObject Reductions.")
    stop("Not Support Reduction: ", reduction)
  }
  cell.umap <- object@reductions[[reduction]]@cell.embeddings # DimReduc

  ## runSubcluster #####
  cell.subcluster <- runSubcluster(object = object, n.core = n.core, cell.range = cell.range, cell.size = cell.size, set.seed = set.seed)
  object@meta.data$cell_subc <- "undefine"
  object@meta.data$cell_subc[match(names(cell.subcluster), rownames(object@meta.data))] <- cell.subcluster
  subc.umap <- aggregate(cell.umap, list(subc.cluster = object@meta.data$cell_subc), mean)
  rownames(subc.umap) <- subc.umap[,1]

  # Subcluster - MeanExpr
  cell.subcluster.table <- table(cell.subcluster)
  total.exp.subcluster <- pbmclapply(names(cell.subcluster.table), function(cell.subcluster.name){ # i = 1
    sub.cells <- names(cell.subcluster)[cell.subcluster %in% cell.subcluster.name]
    # FetchData(object = object, cells = sub.cells, var = rownames(object), slot = "data") # slow
    cluster.mat <- object@assays[[assay]]@data[, sub.cells, drop = F]
    rowMeans(cluster.mat)
  }, mc.cores = n.core) %>% do.call(what = "rbind") %>% t()
  colnames(total.exp.subcluster) <- names(cell.subcluster.table)

  # Construct HematoSubCluster
  hemato.subc <- CreateSubclusterObject(object = total.exp.subcluster, meta.data = subc.umap, sc.meta.data = object@meta.data)
  return(hemato.subc)
}



#'
#' @param meta.data data.frame.\cr Meta data information of sub-cluster.
#'
#' @rdname CreateSubclusterObject
#' @export
#' @method CreateSubclusterObject matrix
#'
CreateSubclusterObject.matrix <- function(object, meta.data, ...){
  subc.cell <- colnames(object)
  if(!all(subc.cell == meta.data[,1])) stop("The column names of Object are not match the first column of meta.data.")
  if(!has.rownames(meta.data)) rownames(meta.data) <- meta.data[,1]
  hemato.subc <- HematoSubCluster(data = object, meta.data = meta.data, misc = list(...))

  ##  runCosineSimilarity #####
  message.time("CosineDist Running")
  scell.ref.path <- system.file("extdata", "HematoMap_scell_ref.rds", package = "HematoMap")
  scell.ref <- readRDS(scell.ref.path)
  # scell.ref <- scell.ref[, traj.data$Cell]
  cosine <- runCosineSimilarity(object = object, scell.ref = scell.ref)
  hemato.subc@dist[["cosine"]] <- cosine
  message.time("CosineDist Complete!")

  return(hemato.subc)
}




#' Run standard Seurat workflow
#'
#' @param object Suerat.\cr Seurat object.
#' @param assay character.\cr Assay name. Default is "RNA".
#' @param force.run logical.\cr Force to run standard Seurat workflow. Default is FALSE.
#'
#' @return Return Seurat object with normalizing, scaling, PCA, clustering.
#' @export
#'
runSeuratWorkflow <- function(object, assay = "RNA", force.run = F){
  if(length(object@commands$NormalizeData.RNA) == 0 || force.run)
    object <- NormalizeData(object, assay = assay, normalization.method = "LogNormalize", scale.factor = 10000, verbose = F)
  if(length(object@commands$FindVariableFeatures.RNA) == 0 || length(VariableFeatures(object)) == 0 || force.run)
    object <- FindVariableFeatures(object, assay = assay, selection.method = "vst", nfeatures = 3000, verbose = F)
  if(length(object@commands$ScaleData.RNA) == 0 || force.run)
    object <- ScaleData(object, assay = assay, features = VariableFeatures(object), verbose = F)
  if(length(object@commands$RunPCA.RNA) == 0 || force.run)
    object <- RunPCA(object, assay = assay, features = VariableFeatures(object), seed.use = 42, verbose = F)
  if(length(object@commands$FindNeighbors.RNA) == 0 || length(object@commands$FindClusters) == 0 || force.run){
    object <- FindNeighbors(object, assay = assay, dims = 1:30, verbose = F)
    object <- FindClusters(object, resolution = 1, verbose = F)
  }
  return(object)
}


#'
#' @param cluster.key character.\cr The key of clustering data in Seurat object.
#'   Default is "seurat_clusters".
#' @param cell.range integer.\cr Range of change of cells number in a sub-cluster
#'   when runing Kmean. Default is 20.
#' @param cell.size integer.\cr Number of cells in a sub-cluster when runing Kmean.
#'   Default is 50.
#' @param set.seed integer.\cr Seed number for Kmean to compute sub-cluster of cell cluster.
#'   Default is 1.
#' @param n.core integer.\cr Number of core for parallel job of running MinMax Kmean cluster.
#'   Default is 1.
#'
#' @rdname runSubcluster
#' @export
#' @method runSubcluster Seurat
#'
runSubcluster.Seurat <- function(object, cluster.key = "seurat_clusters", cell.range = 20, cell.size = 50, set.seed = 1, n.core = 1){
  if(is.null(cluster.key)){
    cell.type.ref <- Idents(object)
  }else{
    cell.type.ref <- setNames(object@meta.data[, cluster.key], rownames(object@meta.data))
  }
  if(is.factor(cell.type.ref)){
    cell.cluster <- levels(cell.type.ref)
  }else{
    cell.cluster <- unique(cell.type.ref)
  }
  if(all(nchar(cell.cluster) < 3)){
    names(cell.cluster) <- paste0("C", cell.cluster)
  }else{
    names(cell.cluster) <- cell.cluster
  }

  message.time("Subcluster Computing")
  .local <- function(i) {
    kmeans.subc <- NULL
    cell.name.subc <- names(cell.cluster)[i]
    if(is.null(cell.name.subc)) cell.name.subc <- cell.cluster[i]
    sub.cells <- names(cell.type.ref)[cell.type.ref %in% cell.cluster[i]]

    if( length(sub.cells) < cell.size + cell.range ){
      sub.cluster <- rep(paste0(cell.name.subc, ":0001"), length(sub.cells) )
      names(sub.cluster) <- sub.cells
      kmeans.subc.re <- sub.cluster

    }else{
      # Re-compute reduction of every cluster
      object.sub <- subset(object, cells = sub.cells)
      object.sub <- FindVariableFeatures(object.sub, selection.method = "vst", nfeatures = 3000, verbose = F)
      object.sub <- ScaleData(object.sub, features = VariableFeatures(object = object.sub), verbose = F)
      object.sub <- RunPCA(object.sub, features = VariableFeatures(object = object.sub), verbose = F)
      cluster.mat <- object.sub@reductions$pca@cell.embeddings[, 1:30]

      # 1. Kmean first
      set.seed(set.seed)
      kmeans.step.1 <- kmeans(cluster.mat, centers = ceiling(nrow(cluster.mat) / cell.size), iter.max = max(10, ceiling(nrow(cluster.mat) / cell.size)) )
      cell.kmeans <- kmeans.step.1$cluster
      stat.1 <- table(cell.kmeans)
      kmeans.subc <- cell.kmeans

      # 2. Loop: Kmean for out-of-range cluster
      cell.undefine <- stat.1[ which(stat.1 < cell.size - cell.range | stat.1 >= cell.size + cell.range ) ]
      while ( length(cell.undefine) > 0 & sum(cell.undefine) > cell.size + cell.range ) {
        cluster.mat.sub <- cluster.mat[names(kmeans.subc)[kmeans.subc %in% names(cell.undefine)], ]
        set.seed(set.seed)
        kmeans.step.2 <- kmeans(cluster.mat.sub, centers = ceiling(nrow(cluster.mat.sub) / cell.size),
                                iter.max = max(10, ceiling(nrow(cluster.mat.sub) / cell.size)) )
        cell.kmeans.2 <- kmeans.step.2$cluster + max(kmeans.subc)
        stat.2 <- table(cell.kmeans.2)

        cell.undefine.new <- stat.2[ which(stat.2 < cell.size - cell.range | stat.2 >= cell.size + cell.range ) ]
        if (length(cell.undefine.new) != length(cell.undefine) ) {
          cell.undefine <- cell.undefine.new
        } else {
          cell.undefine <- NULL
        }
        kmeans.subc[names(cell.kmeans.2)] <- cell.kmeans.2
      }

      # 3. Still out-of-range cluster -> combine to one cluster
      stat.3 <- table(kmeans.subc)
      cell.undefine <- stat.3[ which(stat.3 < cell.size/3) ]
      kmeans.subc[kmeans.subc %in% names(cell.undefine)] <- max(kmeans.subc) + 1

      # Named subcluster
      sub <- table(kmeans.subc)
      kmeans.subc.order <- factor(as.character(kmeans.subc), levels = names(sub)[order(sub, decreasing = T)])
      kmeans.subc.re <- sprintf("%s:%04d", cell.name.subc, as.numeric(kmeans.subc.order))
      names(kmeans.subc.re) <- names(kmeans.subc)

      stat.4 <- table(kmeans.subc.re)
      cell.undefine <- stat.4[ which(stat.4 < 20) ]
      if ( length(cell.undefine) > 0 ) {
        kmeans.subc.re[kmeans.subc.re %in% names(cell.undefine)] <- names(stat.4)[1]
      }
    }
    return(kmeans.subc.re)
  }
  cell.subcluster <- pbmclapply(1:length(cell.cluster), mc.cores = n.core, .local) %>% unlist()
  # sometimes the name of output may be prefixed by "value.", when running only 1 mc.cores.
  if(any(grepl("^value.", names(cell.subcluster)), na.rm = T))
    names(cell.subcluster) <- sub("^value.", "", names(cell.subcluster))
  # if(n.core > 1){
  #   cell.subcluster <- pbmclapply(1:length(cell.cluster), mc.cores = n.core, .local) %>% unlist()
  # }else{
  #   cell.subcluster <- lapply(1:length(cell.cluster), .local) %>% unlist()
  # }
  message.time("Subcluster Complete!")
  return(cell.subcluster)
}





#'
#' @param scell.ref matrix.\cr Reference sub-cluster matrix used for evaluating
#'   the similarity of candidate sub-cluster
#'
#' @return Return a object of \code{\link{CosineDist}}.
#'
#' @importFrom Seurat LogNormalize
#' @importMethodsFrom Matrix t
#' @rdname runCosineSimilarity
#' @method runCosineSimilarity matrix
#'
runCosineSimilarity.matrix <- function(object, scell.ref){
  features.enroll <- intersect(rownames(object), rownames(scell.ref))
  merge.exp.mat <- LogNormalize( cbind(object[features.enroll, ], scell.ref[features.enroll, ] ), verbose = F )
  cosine.sim <- computeCosineSim(merge.exp.mat) # slow
  colnames(cosine.sim) <- rownames(cosine.sim) <- colnames(merge.exp.mat)
  cosine.dist <- computeCosineDist(cosine.sim)
  dimnames(cosine.dist) <- dimnames(cosine.sim)
  CosineDist(sim = cosine.sim, dist = cosine.dist, features.enroll = features.enroll, subc.ref = colnames(scell.ref))
}


#' Compute Cosine Theta and Corresponding Cell Like by Cosine Distance
#'
#' @param hemato.subc \code{\link{HematoSubCluster}}\cr HematoSubCluster object.
#' @param cell.subc character.\cr Names of Candidate sub-cluster
#'
#' @return Return a list containing a vector of cosine theta and a matrix of cell like.
#'
computeThetaLike <- function(hemato.subc, cell.subc = rownames(hemato.subc@meta.data)){
  if(is.null(hemato.subc@dist)) stop("No consine")
  cosine <- hemato.subc@dist[["cosine"]]
  cell.subc.ref <- cosine@subc.ref
  names(cell.subc.ref) <- remove.subc.suffix(cell.subc.ref)
  theta <- sapply(cell.subc, function(x) {
    # 这里包含了NKT的数据 - 如果是用排除
    bc <- as.numeric( rowMeans(cosine@dist[cell.subc.ref, paste0("HSC/MPP:000", 1:3) ]) )
    ab <- as.numeric( mean(cosine@dist[paste0("HSC/MPP:000", 1:3), x ]) )
    ac <- as.numeric( cosine@dist[cell.subc.ref, x ] )

    cos.theta <- (ac ** 2 + bc ** 2 - ab ** 2) / (2 * ac * bc)
    cos.theta[which(cos.theta <= -1)] <- -1
    cos.theta[which(cos.theta >= 1)] <- 1
    names(cos.theta) <- cell.subc.ref
    return(acos( cos.theta ) * 180 / pi )
  })
  like.coef <-  sin( theta / 180 * pi  )

  like.coef[cell.subc.ref[names(cell.subc.ref) %in% c("CD8 Tnaive", "CD8 Teff", "CD8 Tex", "CD8 Tdpe", "CD8 Tmpe", "CD4 Tnaive", "CD4 Tem", "CD4 Treg", "NK", "NK-XCL1", "cycling NK", "cycling NK/T")], ] <- 1
  like.mat <- cosine@sim[cell.subc.ref, cell.subc] * like.coef
  like.mat[which(is.na(like.mat))] <- 1
  return(list(theta = theta, like.mat = like.mat))
}


#' Compute Cell Percentage of Predicted Sub-clusters.
#'
#' @param like.mat matrix.\cr Matrix of like.
#' @param top.subc integer.\cr Number of top like sub-clusters for each candidate
#'   sub-cluster. Default is 5.
#' @param ref.subc character.\cr Names of reference sub-cluster.
#'
#' @return Return cell percentage of sub-clusters.
#'
computePredictSubClustePercentage <- function(like.mat, top.subc = 5, ref.subc){
  total.predict.subc <- c()
  for(i in 1:ncol(like.mat)){
    sub <- like.mat[, i]
    sub <- sub[order(sub, decreasing = T)]
    total.predict.subc <- c(total.predict.subc, head(names(sub), top.subc) )
  }
  sub <- table(total.predict.subc)
  ref.subc.num <- as.numeric(sub[match(ref.subc, names(sub))])
  ref.subc.num[which(is.na(ref.subc.num))] <- 0
  # ref.subc.num.log <- log2(ref.subc.num / top.subc + 1)
  ref.subc.percentage <- ref.subc.num / sum(ref.subc.num) * 100
  return(structure(ref.subc.percentage, names = ref.subc))
}



#' Compute Mean of Like Value for Each Sub-clusters
#'
#' @param hemato.subc \code{\link{HematoSubCluster}}\cr HematoSubCluster object.
#' @param cell.subc character.\cr Names of sub-cluster
#' @param top.subc integer.\cr Number of top like sub-clusters for each candidate
#'   sub-cluster. Default is 5.
#'
#' @return Return mean value of like value of each sub-clusters
#'
computeMeanLike <- function(hemato.subc, cell.subc, top.subc = 5){
  like.mat <- computeThetaLike(hemato.subc = hemato.subc, cell.subc = cell.subc)$like.mat
  total.mean.like <- numeric(nrow(like.mat))
  for(i in 1:nrow(like.mat)){
    sub <- like.mat[i, ]
    total.mean.like[i] <- mean(head(sub[order(sub, decreasing = T)], top.subc))
  }
  return(total.mean.like)
}


#' Compute Lasso Score of BMMC in RNA-seq.
#'
#' @param data matrix.\cr Candidate RNA-seq data in acute leukemia cells.
#' @param normalize logical.\cr Whether to normalize RNA-seq data by normalizeQuantiles.
#'   Default is TRUE.
#'
#' @return Return lasso score of BMMC in RNA-seq.
#'
#' @export
#'
computeLassoScore <- function(data, normalize = T){
  if(!is.matrix(data)) data <- as.matrix(data)
  if(normalize){
    id.common <- intersect(rownames(data), rownames(normal.bm.exp))
    if(length(id.common) == 0) stop("The row names of input data do not match with these of normal.bm.exp")
    data <- normalizeQuantiles(as.matrix(cbind(normal.bm.exp[id.common,,drop = F], data[id.common,,drop = F])))
  }

  # lasso.cell <- c("HSC/MPP","LMPP","CLP","CMP","MDP","GMP","CDP","MEP","pro-Mono","CD14 Mono","CD16 Mono","pre-DC","pDC","cDC1","mo-DC","MKP","MK","pro-Ery1","pro-Ery2","Ery","pre-pro-B","Early pro-B","Late pro-B","pre-B","Immature B","Naive B","Memory B 1","Memory B 2","CD8 Tnaive","CD8 Teff","CD8 Tex","CD8 Tdpe","CD8 Tmpe","CD4 Tnaive","CD4 Tem","CD4 Treg","NK","NK-XCL1")
  mat.score <- matrix(0, nrow = ncol(data), ncol = length(lasso.coef))
  rownames(mat.score) <- colnames(data)
  colnames(mat.score) <- names(lasso.coef)

  for (i in 1:length(lasso.coef)) { # i = 1
    cell.coef <- lasso.coef[[ i ]]
    cell.coef <- cell.coef[names(cell.coef) %in% rownames(data)]
    mat.score[, i] <- colSums(data[names(cell.coef), rownames(mat.score), drop = F] * cell.coef)
  }
  if(normalize){
    ctrl.mat <- mat.score[colnames(normal.bm.exp), ,drop = F]
    mat.score <- mat.score[setdiff(rownames(mat.score), colnames(normal.bm.exp)), ,drop = F]
    for (i in 1:ncol(mat.score)) mat.score[, i] <- mat.score[, i] - mean(ctrl.mat[, i])
  }
  class(mat.score) <- c("HematoLasso", class(mat.score))
  return(mat.score)
}



