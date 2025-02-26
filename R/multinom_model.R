#' Fit a multinomial modelto predict annotation using NMF factors as covariates.
#'
#' @param factors a matrix of factors from NMF
#' @param source_annotations a vector of annotations
#'
#' @return A multinomial model object
#'
#' @import glmnet
#' @importFrom stats predict

fit_multinom_model <- function(factors, source_annotations, alpha=0.5){
  message("Fitting prediction model")
  design <- as.data.frame(cbind(annot=source_annotations, factors))
  print(head(design))

  # mod <-  nnet::multinom(annot ~ ., data = design,
  #                  na.action=na.exclude, maxit=1000)
  mod <- cv.glmnet(x=factors, y=source_annotations, data=design,
                   family = "multinomial", type.multinomial = "grouped",
                   alpha=alpha)
  #p.fit <- predict(mod, predictors=design[grepl("NMF", colnames(design))], type='probs')

  return(mod)
}
