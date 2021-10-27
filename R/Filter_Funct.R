

#' PAD Filtering function
#'
#'A filtering function for mapper that projects $$R$^n$ into $R$
#'
#' @param exp_matrix matrix including the fit residual of the dataset and the healthy state model.
#' @param p integer
#' @param k integer
#'
#' @return A numeric vector including the values produced by the function for each sample in the dataset.
#' @export
#'
#' @examples
#' \dontrun{
#' lp_norm_k_powers(disease_state,2,1)}
#'
lp_norm_k_powers <- function(exp_matrix,p,k){
  if(is.matrix(exp_matrix)){
    return(apply(exp_matrix,2,function(x) (sum(abs(x)^p)^(k/p))))
  }else{
    print("A  matrix object should be provided.")
  }
}

#' SurvPAD Filtering function
#'
#' A filtering function for mapper that projects $$R$^n$ into $R$ that
#'
#' @param exp_matrix matrix including the fit residual of the dataset and the healthy state model.
#' @param p integer
#' @param k integer
#' @param cox_data A matrix with the output of the cox_all_genes function that stores the information of all cox proportional hazard model tests for each gene in the dataset.
#'
#' @return A numeric vector including the values produced by the function for each sample in the dataset.
#' @export
#'
#' @examples
#' \dontrun{
#' lp_norm_k_powers(disease_state,2,1)}
#'
lp_norm_k_powers_surv <- function(exp_matrix,p,k,cox_data){
  if(is.matrix(exp_matrix)){
    cox_vector <- cox_data[rownames(exp_matrix),"z"]
    print(length(cox_vector))
    print(dim(exp_matrix))
    exp_matrix <- exp_matrix * cox_vector
    lp_norm <- apply(exp_matrix,2,function(x) (sum(abs(x)^p)^(k/p)))
    return(lp_norm)
    rownames(exp_matrix)
  }else{
    print("A  matrix object should be provided.")
  }
}
