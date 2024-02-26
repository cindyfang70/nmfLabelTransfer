#' Run NMF on a dataset
#'
#' @param data SingleCellExperiment or SpatialExperiment object
#' @param assay string indicating the assay to run NMF on
#' @param k integer indicating the number of factors for NMF
#' @param seed a random seed
#'
#' @return NMF model object
#'
#' @import SingleCellExperiment
#' @import RcppML
#' @export
run_nmf<- function(data, assay, k, seed=1237){
  A <- assays(data, assay)

  if(missing(k)){
    warning("Number of factors for NMF not specified. Using cross-validation to idenitfy optimal number of clusters")
      k <- find_num_factors(A)
  }


  model <- RcppML::nmf(A, k = k, seed=seed)

  return(model)
}

find_num_factors <- function(A, ranks = c(50,100,200)){
  cv <- singlet::cross_validate_nmf(A, ranks = ranks,
                                    n_replicates = 3,
                                    verbose=3)
  num_factors <- singlet::GetBestRank(cv)

  return(num_factors)
}