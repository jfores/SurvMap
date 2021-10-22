#' Flatten normal tissues
#'
#' Given a matrix which contains the expression values of n healthy tissue samples produces the flattened vector matrix
#' as reported in  Disease-specific genomic analysis: identifying the signature of pathologic biology.
#'
#' @param normal_tiss A normal tissue data gene expression matrix.
#'
#' @return A gene expression matrix containing the flattened versión of the vectors.
#' @export
#'
#' @examples
#' normal_tissue_matrix <- matrix(stats::rnorm(36),nrow=6)
#' flatten_normal_tiss(normal_tissue_matrix)
flatten_normal_tiss <- function(normal_tiss){
  df_out <- normal_tiss
  for(i in 1:ncol(normal_tiss)){
    df_out[,i] <- stats::fitted(stats::lm(normal_tiss[,i] ~ 0 + ., data = data.frame(normal_tiss)[,-i]))
  }
  return(df_out)
}


#' Generate disease component matrix.
#'
#' This functions produces a disease component matrix from the complete dataset and the denoised flattened normal data.
#'
#' @param complete_ds Matrix containing the full dataset.
#' @param normal_space Denoised flattened healthy tissue data.
#'
#' @return disease_component a matrix containing the disease component of the complete dataset.
#' @export
#'
#' @examples
#' full_data <- matrix(stats::rnorm(120),col=20)
#' normal_tissue <- full_data[,11:20]
#' normal_space <- flatten_normal_tiss(normal_tissue)
#' disease_component <- generate_disease_component(full_data,normal_space)
generate_disease_component <- function(complete_ds,normal_space){
  disease_component <- complete_ds
  for(i in 1:ncol(complete_ds)){
    print(i)
    disease_component[,i] <- resid(lm(complete_ds[,i] ~ 0 + ., data = data.frame(normal_space)))
  }
  return(disease_component)
}

