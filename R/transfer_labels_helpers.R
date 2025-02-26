source_nmf_and_model_fitting <- function(source, assay, seed, save_nmf, nmf_path, annotationsName, technicalVarName, alpha,...){
  # 1: run NMF on the source dataset
  source_nmf_mod <- run_nmf(data=source, assay=assay, seed=seed, ...)

  if(save_nmf){
    saveRDS(source_nmf_mod, nmf_path)
  }

  source_factors <- t(source_nmf_mod$h)
  colnames(source_factors) <- paste0("NMF", 1:ncol(source_factors))

  # 2: select the important factors based on correlation
  annots <- colData(source)[[annotationsName]] # compute correlation with domains of interest
  factor_annot_cors <- compute_factor_correlations(source_factors, annots)


  # compute correlation with technial variable
  technicalVar <- colData(source)[[technicalVarName]]
  factor_tech_cors <- compute_factor_correlations(source_factors, technicalVar)

  factors_use_names <- identify_factors_representing_annotations(factor_annot_cors, factor_tech_cors)

  # 3: fit multinomial model on source factors
  #factors_use <- source_factors[,factors_use_names]
  multinom_mod <- fit_multinom_model(source_factors, annots, alpha)

  return(list(source_nmf=source_nmf_mod,
              source_factors=source_factors,
              factor_annot_cors=factor_annot_cors,
              factors_use_names=factors_use_names, multinom=multinom_mod))
}

check_targets_validity <- function(target, assay){
  if (is(target,"SpatialExperiment")){
    if (!(assay %in% assayNames(target))){
      stop(sprintf("Assay %s is not found in the target dataset. %s needs to be available in both the source and target datasets.", assay, assay))
    }

    if(!("gene_name" %in% colnames(rowData(target)))){
      stop("Please provide gene symbols in your target dataset as a column named 'gene_name' in rowData.")
    }

  }
}

check_source_validity <- function(source, assay, annotationsName){
  if(is(source, "SpatialExperiment")){
    if (!(assay %in% assayNames(source))){
      stop(sprintf("Assay %s not found in the source dataset. %s needs to be available in both the source and target datasets.", assay, assay))
    }
    if(!(annotationsName %in% colnames(colData(source)))){
      stop(sprintf("%s is not found in the colData of the source dataset.", annotationsName))
    }
    if(!("gene_name" %in% colnames(rowData(source)))){
      stop("Please provide gene symbols in your source dataset as a column named 'gene_name' in rowData.")
    }

  }
}
