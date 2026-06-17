#' Run NMF on a dataset
#'
#' @param data SingleCellExperiment or SpatialExperiment object
#' @param assay string indicating the assay to run NMF on
#' @param k integer indicating the number of factors for NMF
#' @param seed a random seed
#' @param ... additional arguments passed to singlet::run_nmf or singlet::RunNMF
#'
#' @return NMF model object
#'
#' @import SingleCellExperiment
#' @import RcppML
#' @export
run_nmf<- function(data, assay, k=NULL, seed=1237,...){

  #add in checks for k
  message("Running NMF")

  if(is.null(k)){
    warning("Number of factors for NMF not specified. Using cross-validation to idenitfy optimal number of factors.", immediate. = TRUE)
      #k <- find_num_factors(A)
    set.seed(seed)
    model <- run_rank_determination_nmf(data, assay,...)
  }else{
    A <- as.matrix(assay(data, assay))
    set.seed(seed)
    model <- singlet::run_nmf(A, rank=k,...)
  }

  return(model)
}

find_num_factors <- function(A, ranks = c(50,100,200)){
  if(any(ranks >= ncol(A))){
    stop("ranks must be less than the number of columns in A")
    }
  cv <- singlet::cross_validate_nmf(A, ranks = ranks,
                                    n_replicates = 3,
                                    verbose=3)
  num_factors <- singlet::GetBestRank(cv)

  return(num_factors)
}

#' run_rank_determination_nmf
#'
#' @param data A SingleCellExperiment or SpatialExperiment object
#' @param assay string indicating the assay to run NMF on
#' @param ... additional arguments passed to singlet::RunNMF
#'
#' @return a NMF model object
#' @import singlet SingleCellExperiment
run_rank_determination_nmf <- function(data, assay,...){
  data_nmf <- RunNMF(data, assay=assay,...)
  nmf_mod <- S4Vectors::metadata(data_nmf)$nmf_model
  return(nmf_mod)
}

#' Project source NMF loadings onto a target dataset
#'
#' @param source A SingleCellExperiment or SpatialExperiment used to fit the NMF.
#' @param target A SingleCellExperiment or SpatialExperiment to project.
#' @param assay Name of the assay shared by source and target.
#' @param nmf_model An NMF model with `$w` (gene loadings) and `$d` (scaling).
#' @param harmonize Feature harmonization applied to the shared genes before
#'   projection. `"none"` (default) reproduces the original behaviour. `"zscore"`
#'   standardizes each shared gene and then rescales it to the source's per-gene
#'   mean and standard deviation, so the projection input is placed on the scale
#'   the loadings `w` were trained on while removing platform-specific per-gene
#'   scale/offset differences. This keeps the projected factor scores comparable
#'   to the source factor scores the multinomial model was trained on.
#'
#' @return A cells-by-factors matrix of projected NMF scores.
project_factors <- function(source, target, assay, nmf_model, harmonize = c("none", "zscore")){
  harmonize <- match.arg(harmonize)
  if(is(target, "SpatialExperiment")){
    if (!(assay %in% assayNames(source))){
      stop(sprintf("Assay %s not found in the source dataset. %s needs to be available in both the source and target datasets.", assay, assay))
    }

    if(!("gene_name" %in% colnames(rowData(target)))){
      stop("Please provide gene symbols in your target dataset as a column named 'gene_name' in rowData.")
    }
  }

  if(is(source, "SpatialExperiment")){
    if (!(assay %in% assayNames(source))){
      stop(sprintf("Assay %s not found in the source dataset. %s needs to be available in both the source and target datasets.", assay, assay))
    }

    if(!("gene_name" %in% colnames(rowData(source)))){
      stop("Please provide gene symbols in your source dataset as a column named 'gene_name' in rowData.")
    }
  }



  loadings <- nmf_model$w
  #subset to the genes shared between the source and the target
  i<-intersect(rowData(target)$gene_name, rowData(source)$gene_name) # need to check for gene names in source object later

  if(length(i) == 0){
    stop("No intersecting genes between target and source dataset.")
  }

  rownames(loadings) <- rowData(source)$gene_name
  loadings<-loadings[rownames(loadings) %in% i,]
  loadings <- loadings[unique(rownames(loadings)),] # genes may get duplicated

  target <- target[rowData(target)$gene_name %in% i, ]
  loadings<-loadings[match(rowData(target)$gene_name,rownames(loadings)),]
  #print(any(is.na(loadings)))

  A <- as.matrix(assay(target, assay))

  if(harmonize == "zscore"){
    A <- harmonize_to_source(A, source, assay, gene_order = rowData(target)$gene_name)
  }

  options(RcppML.threads = 0) #line below doesn't work otherwise
  proj<-RcppML::project(data=A, w=loadings, threads=0, L1=0, mask=NULL)

  #print(head(proj))
  factors <- t(proj)/nmf_model$d #scale by the constant factor
  colnames(factors) <- paste0("NMF", 1:ncol(factors))

  return(factors)
}

#' Residualize platform effects across source and target(s) before NMF
#'
#' Restricts source and target(s) to their shared genes (matched by
#' `rowData$gene_name`), stacks them into one matrix, and removes a per-gene
#' platform location effect by fitting `expression ~ platform` for each gene
#' (platform = dataset of origin) and keeping the residuals re-centred on the
#' grand per-gene mean. This is equivalent to `limma::removeBatchEffect` with a
#' single batch factor and no covariates. Negative corrected values are clamped
#' to zero so the result is a valid non-negative input for NMF.
#'
#' Unlike projection-time z-scoring, this couples the source and target(s): NMF
#' must be (re)fit on the corrected source, so the resulting model is specific to
#' the target(s) supplied here and is not a reusable artifact. The single-factor
#' model `~ platform` removes per-gene *location* differences only; per-gene
#' scale/contrast differences are left intact.
#'
#' @param source A SingleCellExperiment/SpatialExperiment used to fit the NMF.
#' @param targets A list of SingleCellExperiment/SpatialExperiment objects.
#' @param assay The assay shared by source and targets.
#'
#' @return A list with `source` (corrected) and `targets` (list of corrected
#'   objects), all subset to the shared genes in a common order.
#' @keywords internal
harmonize_platform_residual <- function(source, targets, assay){
  objs <- c(list(source), targets)

  # shared genes (matched on gene_name) across all datasets
  gene_lists <- lapply(objs, function(o) unique(rowData(o)$gene_name))
  shared <- Reduce(intersect, gene_lists)
  if(length(shared) == 0){
    stop("No genes shared across source and all target datasets; cannot residualize on platform.")
  }

  # align every object to the shared genes in a common order
  align <- function(o){
    idx <- match(shared, rowData(o)$gene_name) # first row per shared gene
    o[idx, ]
  }
  objs <- lapply(objs, align)

  mats <- lapply(objs, function(o) as.matrix(assay(o, assay)))
  n_per <- vapply(mats, ncol, integer(1))
  platform <- factor(rep(paste0("dataset", seq_along(mats)), n_per))

  C <- do.call(cbind, mats)
  grand <- rowMeans(C)
  for(p in levels(platform)){
    idx <- which(platform == p)
    C[, idx] <- C[, idx, drop = FALSE] - rowMeans(C[, idx, drop = FALSE]) + grand
  }
  C[C < 0] <- 0 # keep NMF input non-negative

  # write corrected blocks back into the aligned objects
  ends <- cumsum(n_per); starts <- c(1, head(ends, -1) + 1)
  for(i in seq_along(objs)){
    assay(objs[[i]], assay) <- C[, starts[i]:ends[i], drop = FALSE]
  }

  list(source = objs[[1]], targets = objs[-1])
}

#' Standardize a target expression matrix to the source per-gene distribution
#'
#' Each gene in the target matrix is z-scored using its own mean/sd and then
#' rescaled to the source gene's mean/sd. This removes platform-specific per-gene
#' scale and offset differences while keeping the result on the scale the NMF
#' loadings were learned on. Genes absent from the source, or with zero variance,
#' are passed through (target-only z-score) rather than dropped.
#'
#' @param A A genes-by-cells target expression matrix (rows ordered by `gene_order`).
#' @param source A SingleCellExperiment/SpatialExperiment used to fit the NMF.
#' @param assay The assay name shared by source and target.
#' @param gene_order Character vector of gene symbols giving the row order of `A`.
#'
#' @return A genes-by-cells matrix on the source per-gene scale.
#' @keywords internal
harmonize_to_source <- function(A, source, assay, gene_order){
  src <- as.matrix(assay(source, assay))
  rownames(src) <- rowData(source)$gene_name
  src <- src[!duplicated(rownames(src)), , drop = FALSE]

  src_mu <- rep(0, length(gene_order))
  src_sd <- rep(1, length(gene_order))
  present <- gene_order %in% rownames(src)
  src_mu[present] <- rowMeans(src[gene_order[present], , drop = FALSE])
  src_sd[present] <- apply(src[gene_order[present], , drop = FALSE], 1, stats::sd)

  tgt_mu <- rowMeans(A)
  tgt_sd <- apply(A, 1, stats::sd)

  # guard against zero-variance genes (no scaling information)
  tgt_sd[!is.finite(tgt_sd) | tgt_sd == 0] <- 1
  src_sd[!is.finite(src_sd) | src_sd == 0] <- 1

  (A - tgt_mu) / tgt_sd * src_sd + src_mu
}
