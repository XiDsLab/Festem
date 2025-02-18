
FestemCore <- function(counts,cluster.labels, batch.id,
                            prior.weight = 0.05, prior.weight.filter = 0.9,
                            earlystop = 1e-4, outlier_cutoff = 0.90,
                            min.percent = 0.01,min.cell.num = 30,
                            seed = 321,num.threads = 2,block_size = 40000){
  message("Running Festem ...")
  set.seed(seed)
  cluster.labels <- factor(cluster.labels)
  levels(cluster.labels) <- 1:nlevels(cluster.labels)
  batch.id <- factor(batch.id)
  
  RNGversion("4.4.0")
  cl <- parallel::makeCluster(getOption("cl.cores", num.threads))
  parallel::clusterSetRNGStream(cl, iseed = seed)
  em.result <- vector("list",nlevels(batch.id))
  em.result.f <- vector("list",nlevels(batch.id))
  for (B in 1:nlevels(batch.id)){
    counts.tmp <- counts[,batch.id==(levels(batch.id)[B])]
    cluster.labels.tmp <- cluster.labels[batch.id==(levels(batch.id)[B])]
    # Removing those clusters with no more than 10 cells in this batch
    cluster.labels.tmp <- factor(cluster.labels.tmp)
    counts.tmp <- counts.tmp[,summary(cluster.labels.tmp)[cluster.labels.tmp]>10]
    cluster.labels.tmp <- cluster.labels.tmp[summary(cluster.labels.tmp)[cluster.labels.tmp]>10]
    cluster.labels.tmp <- factor(cluster.labels.tmp)
    levels(cluster.labels.tmp) <- 1:nlevels(cluster.labels.tmp)
    
    ## Outlier
    message(paste0("Batch ",levels(batch.id)[B],": preprocessing ..."))
    if (block_size > nrow(counts.tmp)){
      if (requireNamespace("pbapply", quietly = TRUE)) {
        counts.tmp <- t(pbapply::pbapply(counts.tmp,1,rm.outlier,percent = 0.95, cl = cl))
      } else {
        counts.tmp <- t(parallel::parApply(cl,counts.tmp,1,rm.outlier,percent = 0.95))
      }
    } else{
      block_num <-  ceiling(nrow(counts.tmp) / block_size)
      for (block_iter in 1:block_num){
        if (block_iter == 1){
          if (requireNamespace("pbapply", quietly = TRUE)) {
            counts.tmp.tmp <- pbapply::pbapply(counts.tmp[1:block_size,],1,rm.outlier,percent = 0.95, cl = cl)
          } else {
            counts.tmp.tmp <- parallel::parApply(cl,counts.tmp[1:block_size,],1,rm.outlier,percent = 0.95)
          }
        } else{
          if (requireNamespace("pbapply", quietly = TRUE)) {
            index_end = min(block_size * block_iter, nrow(counts.tmp))
            counts.tmp.tmp <- cbind(counts.tmp.tmp,
                                    pbapply::pbapply(counts.tmp[(block_size * (block_iter-1)+1):index_end,],1,rm.outlier,percent = 0.95, cl = cl))
          } else {
            index_end = min(block_size * block_iter, nrow(counts.tmp))
            counts.tmp.tmp <- cbind(counts.tmp.tmp,
                                    parallel::parApply(cl,counts.tmp[(block_size * (block_iter-1)+1):index_end,],1,rm.outlier,percent = 0.95))
          }
        }
        gc(verbose = FALSE)
      }
      counts.tmp <- t(counts.tmp.tmp)
      rm(counts.tmp.tmp)
    }
    
    rownames(counts.tmp) <- rownames(counts)
    
    ## Sub-sampling
    library.size <- edgeR::calcNormFactors(counts.tmp)
    if (requireNamespace("pbapply", quietly = TRUE)) {
      counts.tmp <- pbapply::pbapply(rbind(library.size,counts.tmp),2,sub.sample, cl = cl)
    } else {
      counts.tmp <- parallel::parApply(cl,rbind(library.size,counts.tmp),2,sub.sample)
    }
    
    if (block_size >= nrow(counts.tmp)){
      if (requireNamespace("pbapply", quietly = TRUE)) {
        counts.tmp <- t(pbapply::pbapply(counts.tmp,1,rm.outlier,percent = 0.95, cl = cl))
      } else {
        counts.tmp <- t(parallel::parApply(cl,counts.tmp,1,rm.outlier,percent = 0.95))
      }
    } else{
      block_num <-  ceiling(nrow(counts.tmp) / block_size)
      for (block_iter in 1:block_num){
        if (block_iter == 1){
          if (requireNamespace("pbapply", quietly = TRUE)) {
            counts.tmp.tmp <- pbapply::pbapply(counts.tmp[1:block_size,],1,rm.outlier,percent = outlier_cutoff, cl = cl)
          } else {
            counts.tmp.tmp <- parallel::parApply(cl,counts.tmp[1:block_size,],1,rm.outlier,percent = outlier_cutoff)
          }
        } else{
          if (requireNamespace("pbapply", quietly = TRUE)) {
            index_end = min(block_size * block_iter, nrow(counts.tmp))
            counts.tmp.tmp <- cbind(counts.tmp.tmp,
                                    pbapply::pbapply(counts.tmp[(block_size * (block_iter-1)+1):index_end,],1,rm.outlier,percent = outlier_cutoff, cl = cl))
          } else {
            index_end = min(block_size * block_iter, nrow(counts.tmp))
            counts.tmp.tmp <- cbind(counts.tmp.tmp,
                                    parallel::parApply(cl,counts.tmp[(block_size * (block_iter-1)+1):index_end,],1,rm.outlier,percent = outlier_cutoff))
          }
        }
        gc(verbose = FALSE)
      }
      counts.tmp <- t(counts.tmp.tmp)
      rm(counts.tmp.tmp)
    }
    # counts.tmp <- t(parallel::parApply(cl,counts.tmp,1,rm.outlier,percent = outlier_cutoff))
    rownames(counts.tmp) <- rownames(counts)
    
    # Only consider genes with non-zero expression in at least min.cell.num
    nonzeros.num <- function(x){sum(x!=0)}
    tmp <- apply(counts.tmp, 1, nonzeros.num)
    min.cell <- min(min.percent*ncol(counts.tmp),min.cell.num)
    #sum(tmp>=min.cell)
    counts.tmp <- counts.tmp[tmp>=min.cell,]
    
    alpha.label <- numeric(nlevels(cluster.labels.tmp)-1)
    for (g in 1:length(alpha.label)) {
      alpha.label[g] <- sum(cluster.labels.tmp==g)/ncol(counts.tmp)
    }
    
    message(paste0("Batch ",levels(batch.id)[B],": EM-test ..."))
    if (requireNamespace("pbapply", quietly = TRUE)) {
      em.result[[B]] <- pbapply::pbapply(counts.tmp,1,em.stat,alpha.ini=rbind(alpha.label,rep(1/nlevels(cluster.labels.tmp),length(alpha.label))),k0=100,C=1e-3,labels = cluster.labels.tmp,
                       group.num = nlevels(cluster.labels.tmp),prior.weight=prior.weight,earlystop = earlystop, cl = cl)
    } else {
      em.result[[B]] <- parallel::parApply(cl,counts.tmp,1,em.stat,alpha.ini=rbind(alpha.label,rep(1/nlevels(cluster.labels.tmp),length(alpha.label))),k0=100,C=1e-3,labels = cluster.labels.tmp,
                                           group.num = nlevels(cluster.labels.tmp),prior.weight=prior.weight,earlystop = earlystop)
    }
    
    message(paste0("Batch ",levels(batch.id)[B],": filtering ..."))
    if (requireNamespace("pbapply", quietly = TRUE)) {
      em.result.f[[B]] <- pbapply::pbapply(counts.tmp,1,em.stat,alpha.ini=rbind(alpha.label,rep(1/nlevels(cluster.labels.tmp),length(alpha.label))),k0=100,C=1e-3,labels = cluster.labels.tmp,
                       group.num = nlevels(cluster.labels.tmp),prior.weight=prior.weight.filter,earlystop = earlystop, cl = cl)
    } else {
      em.result.f[[B]] <- parallel::parApply(cl,counts.tmp,1,em.stat,alpha.ini=rbind(alpha.label,rep(1/nlevels(cluster.labels.tmp),length(alpha.label))),k0=100,C=1e-3,labels = cluster.labels.tmp,
                                             group.num = nlevels(cluster.labels.tmp),prior.weight=prior.weight.filter,earlystop = earlystop)
    }
    
  }
  
  ## EM
  p.tmp <- matrix(NA,nrow = nrow(counts),ncol = length(em.result))
  rownames(p.tmp) <- rownames(counts)
  for (i in 1:length(em.result)){
    p.tmp[colnames(em.result[[i]]),i] <- em.result[[i]][1,]
  }
  p.tmp <- parallel::parApply(cl,p.tmp,1,my.min)
  p.tmp <- p.tmp*length(em.result)
  
  em.tmp <- matrix(NA,nrow = nrow(counts),ncol = length(em.result))
  rownames(em.tmp) <- rownames(counts)
  for (i in 1:length(em.result)){
    em.tmp[colnames(em.result.f[[i]]),i] <- em.result.f[[i]][2,]
  }
  
  em.tmp <- apply(em.tmp,1,my.max)
  parallel::stopCluster(cl)
  
  result.EM = data.frame(names = rownames(counts), p = p.adjust(p.tmp,"BH"), EM = em.tmp)
  tmp.na <- result.EM[is.na(result.EM$p),1]
  result.EM <- result.EM[!is.na(result.EM$p),]
  gene.names <- result.EM[result.EM$p<0.05 & result.EM$EM>0,]
  gene.names <- gene.names[order(-gene.names$EM),]
  
  tmp <- result.EM[result.EM$p>=0.05 & result.EM$EM>0,]
  tmp <- tmp[order(-tmp$EM),]
  gene.names <- rbind(gene.names,tmp)
  tmp <- result.EM[result.EM$EM<=0,]
  tmp <- tmp[order(tmp$p,-tmp$EM),]
  gene.names <- rbind(gene.names,tmp)
  gene.names <- gene.names[,1]
  gene.names <- c(gene.names,tmp.na)
  EM <- gene.names
  
  list(stat = result.EM,genelist = EM)
}

#' Festem feature selection and differential expression gene analysis
#' 
#' @description Run Festem algorithm with Seurat pipeline or original count matrix.
#' 
#' @param object A Seurat object, matrix or dgCMatrix containing the original count matrix. 
#' If a Seurat object is provided, further information on cells (such as pre-clustering and 
#' batches) can be included in the metadata slot and specified with parameters \code{prior} 
#' and \code{batch}. If a matrix or dgCMatrix is provided, its rows should represent different
#' genes and columns stand for cells. We recommend setting the rownames as the name of genes.
#' If row names are not specified, numerical indices will be created to represent genes
#' according to their order in the input matrix.
#' @param ... other parameters
#' @param G An integer specifying number of mixing component (roughly speaking, number of clusters) 
#' used in EM-test. When \code{object} is a Seurat object and \code{G} is null, \code{G} will be automatically determined
#' by running Louvain clustering using parameters specified by \code{prior_parameters}.
#' Empirically, we recommend setting it slightly larger than the expected number
#' of clusters to gain more statistical power. Theoretically speaking, if \code{G} is larger than 
#' the actual number of clusters, the test is still valid and does not lose any power; however, if it 
#' is smaller than the actual number, the power might decrease slightly but the test is still valid 
#' (i.e. type I error can be controlled).
#' @param prior A string or vector specifying prior (pre-clustering label), can be "HVG" (default), "active.ident", the name of a column  
#' in metadata storing the pre-clustering label for each cell or a vector containing labels for each cell.
#' If "HVG", pre-clustering will be automatically performed with parameters specified by \code{prior_parameters}.
#' If a vector is provided, it should have the same length as the number of cells.
#' @param batch A string or vector specifying batches for cells, can be \code{NULL} (default), the name of a column in metadata
#' storing batch indices, or a vector containing batch indices for cells. If \code{NULL}, all cells are assumed
#' coming from the same batch.
#' @param prior.weight Numerical, specifying prior weight used for EM-test (default = 0.05).
#' @param outlier_cutoff Numerical, ranging from 0 to 1. Outliers are identified only in the largest \eqn{\alpha} proportion counts
#' specified by this parameter.
#' @param seed Numerical, seed of the random number generator.
#' @param num.threads Integer, number of CPU cores to use.
#' @param FDR_level Numerical, nomical FDR level to determine number of selected features (i.e. differential expression genes).
#' @param prior.weight.filter Numerical, specifying prior weight used to filter out genes that have expression patterns significantly 
#' conflict with the pre-clustering labels (default = 0.9).
#' @param min.percent,min.cell.num Numerical, specifying the filtering threshold for a gene to be tested. If a gene has non-zero expression
#' in at most \code{min(n*min.percent,min.cell.num)} cells, where \code{n} is the total cell numbers, then this gene will not be tested.
#' @param earlystop Numerical, specifying threshold of absolute changes in penalized log-likelihood (See \sQuote{Details} in \code{\link{em_test}} ) for early-stopping. 
#' Empirically, if the data is far from homogeneously distributed, the penalized log-likelihood will increase rapidly in the early stage of 
#' iterations. Later iterations only increase the penalized likelihood marginally. Therefore, specifying an early-stopping threshold can save 
#' some time while still guarantee the validity of p-values.
#' @param block_size Integer, number of genes processed at a time when removing outliers. If parApply returns an out-of-memory error, try use a smaller block size.
#' 
#' @examples 
#' # Load example data
#' data("example_data")
#' example_data <- Seurat::CreateSeuratObject(example_data$counts,meta.data = example_data$metadata)
#' example_data <- RunFestem(example_data,2,prior = "labels")
#' Seurat::VariableFeatures(example_data)[1:10]
#' # Use smaller block size for memory issue
#' data("example_data")
#' example_data <- Seurat::CreateSeuratObject(example_data$counts,meta.data = example_data$metadata)
#' example_data <- RunFestem(example_data,2,prior = "labels",block_size = 50)
#' Seurat::VariableFeatures(example_data)[1:10]
#' @details For details of Festem, see Chen, Wang, et al (2023).
#' @seealso [em_test()]
#' @references Chen, Z., Wang, C., Huang, S., Shi, Y., & Xi, R. (2024). Directly selecting cell-type marker genes for single-cell clustering analyses. Cell Reports Methods, 4(7).
#' @rdname RunFestem
#' @export


RunFestem <- function(object,...){
  UseMethod("RunFestem")
}

#' @rdname RunFestem
#' @param prior_parameters A list containing parameters to do pre-clustering with HVG. \code{HVG_num} specified
#' number of HVG to use and \code{PC_dims} is the number of PC dimensions to use.
#' @param assay Name of Assay used
#' @return Seurat (version 4) object. Selected features placed into
#' \code{var.features} object. Test statistics, p-values and ranks for genes are placed into columns of
#' \code{meta.data} in the slot, with column names "p", "EM", and "Festem_rank", respectively.  For downstream Seurat analyses,
#' use default variable features.
#' @export
RunFestem.Seurat <- function(object,G = NULL,prior = "HVG", batch = NULL,
                             prior.weight = 0.05, prior.weight.filter = 0.9,
                             earlystop = 1e-4, outlier_cutoff = 0.90,
                             min.percent = 0.01,min.cell.num = 30,
                             seed = 321,num.threads = 1,FDR_level = 0.05,block_size = 40000, assay = "RNA",
                             prior_parameters = list(HVG_num = 2000,PC_dims = 50, resolution = 0.7),...){
  if (!requireNamespace('Seurat', quietly = TRUE)) {
    stop("Running Festem on a Seurat object requires Seurat")
  }
  object@active.assay <- assay
  # prior: "HVG", "active.ident", column name of meta.data, vector
  # batch: NULL, column name of meta.data, vector
  # Otherwise, prior will be generated with top 2000 HVGs and all cells will be taken as one single batch.
  if (prior == "HVG"){
    message("Generating priors ... ")
    object <- Seurat::FindVariableFeatures(object = object, selection.method = "vst", 
                                   nfeatures = prior_parameters[["HVG_num"]] %or% 8000,
                                   verbose = FALSE)
    object <- Seurat::NormalizeData(object)
    object <- Seurat::ScaleData(object)
    object <- Seurat::RunPCA(object,verbose = FALSE)
    object <- Seurat::FindNeighbors(object = object, dims = 1:(prior_parameters[["PC_dims"]] %or% 50))
    if (is.null(G)){
      object <- Seurat::FindClusters(object, resolution = prior_parameters[["resolution"]], verbose = FALSE)
      G <- nlevels(object@active.ident)
      prior_label <- object@active.ident
    } else{
      for (k in 1:100){
        object <- Seurat::FindClusters(object, resolution = 0.01*k, verbose = FALSE)
        if (nlevels(object@active.ident)>=G){
          if (nlevels(object@active.ident)==G){
            prior_label <- object@active.ident
            break
          } else{
            stop(paste0("Cannot generate a prior with ",G," clusters. Try setting G as ",nlevels(object@active.ident)," instead and re-run Festem."))
          }
        }
      }
    }
  } else if (prior == "active.ident"){
    prior_label <- object@active.ident
  } else if (length(prior) == 1){
    if (prior %in% colnames(object@meta.data)){
      prior_label <- object@meta.data[,prior]
      if (length(unique(prior_label)) != G){
        warning("Number of clusters in the prior do not equal to G. We take number of clusters in the prior as the number of mixing components in the EM-test. If G is the number of mixing components wanted, please use another prior with number of clusters equals to G.")
        G <- length(unique(prior_label))
      }
    } else{
      stop(paste0(prior," was not found in metadata"))
    }
  } else if (length(prior) == ncol(object)){
    prior_label <- factor(prior)
    if (nlevels(prior_label) != G){
      warning("Number of clusters in the prior do not equal to G. We take number of clusters in the prior as the number of mixing components in the EM-test. If G is the number of mixing components wanted, please use another prior with number of clusters equals to G.")
      G <- nlevels(prior_label)
    }
  } else{
    stop("Prior must be 'HVG', 'active.ident', a column name of metadata or a vector of labels for cells. If a vector of labels were provided, please check whether its length is the same as the number of cells.")
  }
  
  if (is.null(batch)){
    batch_id <- rep(1,ncol(object))
  } else if (length(batch) == 1){
    if (batch %in% colnames(object@meta.data)){
      batch_id <- object@meta.data[,batch]
    } else{
      stop(paste0(batch," was not found in metadata"))
    }
  } else if (length(batch) == ncol(object)){
    batch_id <- batch
  } else{
    stop("Batch must be NULL, a column name of metadata or a vector of batch indices for cells. If a vector were provided, please check whether its length is the same as the number of cells.")
  }
  
  if (G==1){
    stop("Number of components must be larger than 1.")
  }
  
  counts_use <- SeuratObject::LayerData(object,assay = assay,layer = "counts")
  
  ## Filter out batches with too few data, since we require that in each batch there are at least two clusters with more than 10 cells
  filter_table <- table(prior_label,batch_id)
  filter_criteria <- apply(filter_table,2,function(x){y <- sort(x,decreasing = T);y[2]})
  batch_keep <- names(filter_criteria)[filter_criteria > 10]
  filter_flag <- batch_id %in% batch_keep
  counts_use <- counts_use[,filter_flag]
  prior_label <- prior_label[filter_flag]
  batch_id <- batch_id[filter_flag]
  
  print(paste0("The following batches are removed due to small sample size: ",
               paste0(names(filter_criteria)[filter_criteria < 10],collapse = ", ")))
  
  result <- FestemCore(counts_use,cluster.labels = prior_label, batch.id = batch_id,
                             prior.weight = prior.weight, prior.weight.filter = prior.weight.filter,
                             earlystop = earlystop, outlier_cutoff = outlier_cutoff,
                             min.percent = min.percent,min.cell.num = min.cell.num,
                             seed = seed,num.threads = num.threads,block_size = block_size)
  num_features <- sum(result[["stat"]]$p<FDR_level & result[["stat"]]$EM>0)
  Seurat::VariableFeatures(object) <- result[["genelist"]][1:min(num_features,6000)]
  
  meta_to_add <- data.frame("Festem_p_adj" = rep(NA,nrow(object[[assay]])),
                            "Festem_EM" = rep(NA,nrow(object[[assay]])),
                            "Festem_rank" = rep(NA,nrow(object[[assay]])))
  rownames(meta_to_add) <- rownames(object)
  meta_to_add[result[["stat"]]$names,"Festem_p_adj"] <- result[["stat"]]$p
  meta_to_add[result[["stat"]]$names,"Festem_EM"] <- result[["stat"]]$EM
  meta_to_add[,"Festem_rank"] <- sapply(rownames(object),function(x){
    tmp <- which(x==result[["genelist"]])
    if (length(tmp)>0){
      tmp
    } else{
      NA
    }
  })
  object[[assay]] <- SeuratObject::AddMetaData(object[[assay]],meta_to_add)
                                                          
  return(object)
}

#' @rdname RunFestem
#' @return A list containing test results, gene ranks and selected gene list. The first element of 
#' the list is a data.frame whose columns are gene names, p-values and EM statistics for filtering.
#' The second element is a vector of genes ordered according to EM-test. The third element is a 
#' vector of selected genes.
#' @export
RunFestem.matrix <- function(object,G,prior = NULL, batch = NULL,
                             prior.weight = 0.05, prior.weight.filter = 0.9,
                             earlystop = 1e-4, outlier_cutoff = 0.90,
                             min.percent = 0.01,min.cell.num = 30,block_size = 40000,
                             seed = 321,num.threads = 1,FDR_level = 0.05,...){
  if (is.null(rownames(object))){
    rownames(object) <- paste0("Gene ",1:nrow(object))
  }
  
  if (prior.weight > 0 & is.null(prior)){
    stop("Please provide a valid prior such as pre-clustering labels.")
  }
  
  if (prior.weight == 0){
    prior_label <- sample(1:G,ncol(object),replace = T)
  } else{
    if (length(prior) == ncol(object)){
      prior_label <- factor(prior)
      if (nlevels(prior_label) != G){
        warning("Number of clusters in the prior do not equal to G. We take number of clusters in the prior as the number of mixing components in the EM-test. If G is the number of mixing components wanted, please use another prior with number of clusters equals to G.")
        G <- nlevels(prior_label)
      }
    } else{
      stop("Prior must be a vector of labels for cells, please check whether its length is the same as the number of cells.")
    }
  }
  
  if (is.null(batch)){
    batch_id <- rep(1,ncol(object))
  }  else if (length(batch) == ncol(object)){
    batch_id <- batch
  } else{
    stop("Batch must be NULL, a column name of metadata or a vector of batch indices for cells. If a vector were provided, please check whether its length is the same as the number of cells.")
  }
  
  ## Filter out batches with too few data, since we require that in each batch there are at least two clusters with more than 10 cells
  filter_table <- table(prior_label,batch_id)
  filter_criteria <- apply(filter_table,2,function(x){y <- sort(x,decreasing = T);y[2]})
  batch_keep <- names(filter_criteria)[filter_criteria > 10]
  filter_flag <- batch_id %in% batch_keep
  object <- object[,filter_flag]
  prior_label <- prior_label[filter_flag]
  batch_id <- batch_id[filter_flag]
  
  print(paste0("The following batches are removed due to small sample size: ",
               paste0(names(filter_criteria)[filter_criteria < 10],collapse = ", ")))
  
  result <- FestemCore(object,cluster.labels = prior_label, batch.id = batch_id,
                             prior.weight = prior.weight, prior.weight.filter = prior.weight.filter,
                             earlystop = earlystop, outlier_cutoff = outlier_cutoff,
                             min.percent = min.percent,min.cell.num = min.cell.num,
                             seed = seed,num.threads = num.threads,block_size = block_size)
  num_features <- sum(result[["stat"]]$p<FDR_level & result[["stat"]]$EM>0)
  result[["clustering_features"]] <- result[["genelist"]][1:min(num_features,6000)]
  names(result)[2] <- "generank"
  return(result)
}

#' @rdname RunFestem
#' @return A list containing test results, gene ranks and selected gene list. The first element of 
#' the list is a data.frame whose columns are gene names, p-values and EM statistics for filtering.
#' The second element is a vector of genes ordered according to EM-test. The third element is a 
#' vector of selected genes.
#' @export
RunFestem.Matrix <- function(object,G,prior = NULL, batch = NULL,
                             prior.weight = 0.05, prior.weight.filter = 0.9,
                             earlystop = 1e-4, outlier_cutoff = 0.90,
                             min.percent = 0.01,min.cell.num = 30,block_size = 40000,
                             seed = 321,num.threads = 1,FDR_level = 0.05,...){
  if (is.null(rownames(object))){
    rownames(object) <- paste0("Gene ",1:nrow(object))
  }
  
  if (prior.weight > 0 & is.null(prior)){
    stop("Please provide a valid prior such as pre-clustering labels.")
  }
  
  if (prior.weight == 0){
    prior_label <- sample(1:G,ncol(object),replace = T)
  } else{
    if (length(prior) == ncol(object)){
      prior_label <- factor(prior)
      if (nlevels(prior_label) != G){
        warning("Number of clusters in the prior do not equal to G. We take number of clusters in the prior as the number of mixing components in the EM-test. If G is the number of mixing components wanted, please use another prior with number of clusters equals to G.")
        G <- nlevels(prior_label)
      }
    } else{
      stop("Prior must be a vector of labels for cells, please check whether its length is the same as the number of cells.")
    }
  }
  
  if (is.null(batch)){
    batch_id <- rep(1,ncol(object))
  }  else if (length(batch) == ncol(object)){
    batch_id <- batch
  } else{
    stop("Batch must be NULL, a column name of metadata or a vector of batch indices for cells. If a vector were provided, please check whether its length is the same as the number of cells.")
  }
  
  ## Filter out batches with too few data, since we require that in each batch there are at least two clusters with more than 10 cells
  filter_table <- table(prior_label,batch_id)
  filter_criteria <- apply(filter_table,2,function(x){y <- sort(x,decreasing = T);y[2]})
  batch_keep <- names(filter_criteria)[filter_criteria > 10]
  filter_flag <- batch_id %in% batch_keep
  object <- object[,filter_flag]
  prior_label <- prior_label[filter_flag]
  batch_id <- batch_id[filter_flag]
  
  print(paste0("The following batches are removed due to small sample size: ",
               paste0(names(filter_criteria)[filter_criteria < 10],collapse = ", ")))
  
  result <- FestemCore(object,cluster.labels = prior_label, batch.id = batch_id,
                       prior.weight = prior.weight, prior.weight.filter = prior.weight.filter,
                       earlystop = earlystop, outlier_cutoff = outlier_cutoff,
                       min.percent = min.percent,min.cell.num = min.cell.num,
                       seed = seed,num.threads = num.threads,block_size = block_size)
  num_features <- sum(result[["stat"]]$p<FDR_level & result[["stat"]]$EM>0)
  result[["clustering_features"]] <- result[["genelist"]][1:min(num_features,6000)]
  names(result)[2] <- "generank"
  return(result)
}