#' Fit a multinomial model to predict annotation using NMF factors (and an
#' optional covariate) as predictors.
#'
#' @param factors a matrix of factors from NMF (cells x factors).
#' @param source_annotations a vector of annotations (the response).
#' @param alpha elasticnet mixing parameter, with 0 <= alpha <= 1.
#' @param covariate Optional vector of an additional covariate (e.g. diagnosis),
#'   one value per cell, aligned to the rows of `factors`. If supplied it is
#'   added to the design matrix alongside the NMF factors. A factor/character
#'   covariate is expanded into dummy columns; a numeric covariate is used as-is.
#' @param penalize_covariate If FALSE (default), the covariate columns are left
#'   unpenalized (`penalty.factor = 0`) so the model always adjusts for them.
#'
#' @return A `cv.glmnet` model object. When a covariate is used, the encoding is
#'   stored on the model via `attr(mod, "covariate_meta")` so the identical
#'   design can be reconstructed for the target at prediction time;
#'   `attr(mod, "n_factors")` records how many leading columns are NMF factors.
#'
#' @import glmnet
#' @importFrom stats predict model.matrix
fit_multinom_model <- function(factors, source_annotations, alpha=0.5,
                               covariate=NULL, penalize_covariate=FALSE){
  message("Fitting prediction model")

  x <- as.matrix(factors)
  n_factors <- ncol(x)
  cov_meta <- NULL

  if(!is.null(covariate)){
    if(anyNA(covariate)){
      stop("`covariate` contains NAs; remove or impute them before fitting.")
    }
    cov_meta <- covariate_meta(covariate)
    x <- cbind(x, build_covariate_design(covariate, cov_meta))
    message(sprintf("Including covariate in the model (%d column(s))", ncol(x) - n_factors))
  }

  penalty <- rep(1, ncol(x))
  if(!is.null(covariate) && !penalize_covariate){
    penalty[(n_factors + 1):ncol(x)] <- 0   # always retain the covariate
  }

  mod <- cv.glmnet(x = x, y = source_annotations,
                   family = "multinomial", type.multinomial = "grouped",
                   alpha = alpha, penalty.factor = penalty)

  attr(mod, "n_factors") <- n_factors
  attr(mod, "covariate_meta") <- cov_meta
  return(mod)
}

#' Capture the encoding of a covariate so it can be reproduced on new data
#'
#' @param covariate A vector (numeric, factor, or character).
#' @return A list describing how to build the design matrix for this covariate.
#' @keywords internal
covariate_meta <- function(covariate){
  if(is.numeric(covariate)){
    return(list(type = "numeric"))
  }
  lev <- levels(factor(covariate))
  if(length(lev) < 2){
    stop("A factor/character covariate must have at least 2 levels.")
  }
  list(type = "factor", levels = lev)
}

#' Build a covariate design matrix consistent with a stored encoding
#'
#' Used both when fitting (to set the encoding) and when predicting (to
#' reproduce the exact same columns for the target).
#'
#' @param covariate A vector of covariate values.
#' @param meta The encoding returned by [covariate_meta].
#' @return A numeric matrix with one row per element of `covariate`.
#' @keywords internal
build_covariate_design <- function(covariate, meta){
  if(identical(meta$type, "numeric")){
    m <- matrix(as.numeric(covariate), ncol = 1)
    colnames(m) <- "covariate"
    return(m)
  }
  xf <- factor(covariate, levels = meta$levels)
  if(anyNA(xf)){
    unseen <- setdiff(unique(as.character(covariate)), meta$levels)
    stop("Covariate contains levels not present in the training data: ",
         paste(unseen, collapse = ", "))
  }
  mm <- stats::model.matrix(~xf)[, -1, drop = FALSE]
  colnames(mm) <- paste0("covariate", meta$levels[-1])
  mm
}

#' Append the target covariate columns to a projection matrix before prediction
#'
#' If the model was fit without a covariate, the projections are returned
#' unchanged. Otherwise the covariate is read from `colData(target)[[name]]`,
#' encoded the same way as during training, and column-bound to the projections.
#'
#' @param projections A cells-by-factors matrix of projected NMF scores.
#' @param mod A model from [fit_multinom_model].
#' @param target The target SingleCellExperiment/SpatialExperiment.
#' @param diagnosisName Name of the covariate column in `colData(target)`.
#' @return A matrix suitable for `predict(mod, newx = ...)`.
#' @keywords internal
augment_with_covariate <- function(projections, mod, target, diagnosisName){
  meta <- attr(mod, "covariate_meta")
  if(is.null(meta)){
    return(projections)
  }
  if(is.null(diagnosisName) || !(diagnosisName %in% colnames(colData(target)))){
    stop("The model was trained with a covariate, so the target must contain a ",
         "colData column named '", diagnosisName, "'.")
  }
  covariate <- colData(target)[[diagnosisName]]
  if(anyNA(covariate)){
    stop("Target covariate '", diagnosisName, "' contains NAs.")
  }
  cbind(projections, build_covariate_design(covariate, meta))
}
