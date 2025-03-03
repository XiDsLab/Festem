marker_allocate_core <- function(norm.data, type, sig.level = 0.05, plot_result = F, dispersion = NULL){
  # This function returns a logical vector with the same length 
  # of "type" in which each element denotes whether this gene is 
  # highly expressed in this type
  # norm.data should be normal distributed
  # type is a factor denoting which cell belongs to which type
  # "A" means the highest level
  lm.data <- lm(norm.data~0+type)
  suppressWarnings(SK.result <- ScottKnott::SK(lm.data,which = "type",sig.level = sig.level))
  if (plot_result){
    plot(SK.result,dispersion = dispersion)
  }
  clus <- SK.result$out$Result
  clus[,1] <- as.numeric(clus[,1])
  clus <- clus[order(clus[,1],decreasing = T),]
  clus.name <- rownames(clus)
  clus <- clus[,-1]
  if (!is.null(ncol(clus))) clus <- apply(clus, 1, paste0,collapse = "")
  clus <- apply(matrix(clus,nrow = 1), 2, toupper)
  names(clus) <- clus.name
  if (length(clus)>0){
    clus[levels(type)]
  } else{
    clus <- rep("A",nlevels(type))
    names(clus) <- levels(type)
    clus
  }
}

#' Assigns DEGs to the clusters as their markers
#' 
#' @description Given a set of markers (typically provided by \code{\link{RunFestem}}), 
#' use Scott-Knott test to decide in which groups they are up-regulated or down-regulated.
#' 
#' @param object A Seurat object or a matrix or dgCMatrix containing the normalized count matrix. 
#' @param marker A vector containing names of marker genes, for example, output gene list from \code{\link{RunFestem}}.
#' @param ... other parameters
#' @param num_cores Integer, number of CPU cores to use (default = 1).
#' @param FDR_level Numerical, nominal FDR level for multiple testing (default = 0.05).
#' 
#' @returns A gene by cluster matrix indicating expression levels of each gene in different clusters. "A" stands for the highest expression level. We recommend using genes with expression level "A" and "B" as markers.
#' 
#' @examples 
#' # Load example data
#' data("example_data")
#' example_data <- Seurat::CreateSeuratObject(example_data$counts,meta.data = example_data$metadata)
#' example_data <- RunFestem(example_data,2,prior = "labels")
#' result <- AllocateMarker(example_data,rownames(example_data)[1:10],group_by = "labels")
#' print(result)
#' @seealso [RunFestem()]
#' @references A. J. Scott, M. Knott, A cluster analysis method for grouping means in the analysis of variance. Biometrics 30, 507-512 (1974).
#' @rdname AllocateMarker
#' @export

AllocateMarker <- function(object,marker,...){
  UseMethod("AllocateMarker")
}

#' @rdname AllocateMarker
#' 
#' @param group_by A character string or \code{NULL} (default). If \code{NULL}, \code{active.ident} is used
#' If a string, then the corresponding column in \code{meta.data} is taken as labels for cells.
#' @param assay Name of Assay used
#' 
#' @export
AllocateMarker.Seurat <- function(object,marker,group_by = NULL,assay = "RNA", num_cores = 1,
                                  FDR_level = 0.05,debug = F,...){
  if (!requireNamespace('Seurat', quietly = TRUE)) {
    stop("Running Festem on a Seurat object requires Seurat")
  }
  object@active.assay <- assay
  cl <- parallel::makeCluster(getOption("cl.cores", num_cores))
  object <- Seurat::NormalizeData(object)
  if (is.null(group_by)){
    cluster <- object@active.ident
  } else if (group_by %in% colnames(object@meta.data)){
    cluster <- object@meta.data[,group_by]
  } else{
    stop("group_by should be the name of a column in metadata.")
  }
  
  if (length(setdiff(marker,rownames(object)))>0){
    warning(paste0("The following markers are not found: ",
                   paste(setdiff(marker,rownames(object)),collapse = ", "),"."))
    marker <- intersect(marker,rownames(object))
    if (length(marker) == 0){
      stop("At least one valid marker is needed.")
    }
  }
  
  data_use <- SeuratObject::LayerData(object,assay = assay,layer = "data",features = marker)
  if (debug){
      gene.allocation <- apply(data_use, 1, 
                               marker_allocate_core, type = factor(cluster),
                               sig.level = FDR_level)
  }else{
    if (requireNamespace("pbapply", quietly = TRUE)) {
      gene.allocation <- pbapply::pbapply(data_use, 1, 
                                          marker_allocate_core, type = factor(cluster),
                                          sig.level = FDR_level, cl = cl)
    } else {
      gene.allocation <- parallel::parApply(cl,data_use, 1, 
                                            marker_allocate_core, type = factor(cluster),
                                            sig.level = FDR_level)
    }
  }
  parallel::stopCluster(cl)
  
  ### Format output
  message("Formatting output ...")
  gene.allocation <- t(gene.allocation)
  output <- data.frame()
  gene.allocation <- gene.allocation[rowSums(gene.allocation!="A")>0,]
  for (i in 1:ncol(gene.allocation)){
    marker.tmp <- rownames(gene.allocation)[gene.allocation[,i] %in% c("A","B")]
    if (length(marker.tmp) == 0){
      next
    }
    if (is.null(group_by)){
      FC.tmp <- Seurat::FoldChange(object,ident.1 = colnames(gene.allocation)[i],features = marker.tmp)
    } else if (group_by %in% colnames(object@meta.data)){
      FC.tmp <- Seurat::FoldChange(object,ident.1 = colnames(gene.allocation)[i],group.by = group_by,features = marker.tmp)
    } else{
      stop("group_by should be the name of a column in metadata.")
    }
    
    output.tmp <- cbind(FC.tmp,gene = marker.tmp,
                        p_adj = object[[assay]][[]][marker.tmp,"Festem_p_adj"],
                        cluster = colnames(gene.allocation)[i])
    output <- rbind(output,output.tmp)
    output <- output[order(output$avg_log2FC,decreasing = T),]
  }
  
  return(output)
}

#' @rdname AllocateMarker
#' 
#' @param label A vector containing labels for cells.
#' 
#' @export
AllocateMarker.matrix <- function(object,marker,label,num_cores = 1,
                                  FDR_level = 0.05,...){
  # normalized counts are needed
  cl <- parallel::makeCluster(getOption("cl.cores", num_cores))
  if (length(label)!=ncol(object)){
    stop("Labels should have the same length as cell numbers.")
  }
  
  if (length(setdiff(marker,rownames(object)))>0){
    warning(paste0("The following markers are not found: ",
                   paste(setdiff(marker,rownames(object)),collapse = ", "),"."))
    marker <- intersect(marker,rownames(object))
    if (length(marker) == 0){
      stop("At least one valid marker is needed.")
    }
  }
  
  if (requireNamespace("pbapply", quietly = TRUE)) {
    gene.allocation <- pbapply::pbapply(object, 1, 
                                        marker_allocate_core, type = factor(label),
                                        sig.level = FDR_level, cl = cl)
  } else {
    gene.allocation <- parallel::parApply(cl,object, 1, 
                                          marker_allocate_core, type = factor(label),
                                          sig.level = FDR_level)
  }
  parallel::stopCluster(cl)
  return(t(gene.allocation))
}

#' @rdname AllocateMarker
#' 
#' @param label A vector containing labels for cells.
#' 
#' @export
AllocateMarker.Matrix <- function(object,marker,label,num_cores = 1,
                                  FDR_level = 0.05,...){
  # normalized counts are needed
  cl <- parallel::makeCluster(getOption("cl.cores", num_cores))
  if (length(label)!=ncol(object)){
    stop("Labels should have the same length as cell numbers.")
  }
  
  if (length(setdiff(marker,rownames(object)))>0){
    warning(paste0("The following markers are not found: ",
                   paste(setdiff(marker,rownames(object)),collapse = ", "),"."))
    marker <- intersect(marker,rownames(object))
    if (length(marker) == 0){
      stop("At least one valid marker is needed.")
    }
  }
  
  if (requireNamespace("pbapply", quietly = TRUE)) {
    gene.allocation <- pbapply::pbapply(object, 1, 
                                        marker_allocate_core, type = factor(label),
                                        sig.level = FDR_level, cl = cl)
  } else {
    gene.allocation <- parallel::parApply(cl,object, 1, 
                                          marker_allocate_core, type = factor(label),
                                          sig.level = FDR_level)
  }
  parallel::stopCluster(cl)
  return(t(gene.allocation))
}
