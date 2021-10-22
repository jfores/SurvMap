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
