#' Gene selection based on variability
#'
#' @param disease_component_tumors A matrix with the disease component data of the entire dataset.
#' @param percent Percentile for gene selection.
#'
#' @return A vector containing the names of the selected genes.
#' @export
#'
#' @examples
#' \dontrun{
#' gene_selection(disease_component_tumors,0.99)
#' }
gene_selection <- function(disease_component_tumors,percent = 0.85){
  MaxAbs595 <- apply(disease_component_tumors,1,function(x) base::max(base::abs(stats::quantile(x,c(0.05,0.95)))))
  selected_genes <- names(MaxAbs595[MaxAbs595 > stats::quantile(MaxAbs595,percent)])
  return(selected_genes)
}



#' Survival analysis based on gene expression levels.
#'
#' Carries out univariate cox protportional hazard models for the expression levels of each gene included in the dataset and its link with relapse-free or overall survival.
#'
#' @param eData Expression data for disease samples
#' @param time_vector Vector including time to relapse or time to death information
#' @param event_vector Numeric vector indicating if relapse or death have been produced.
#'
#' @return A matrix with the results of the application of proportional hazard models using the expression levels of eahc gene as covariate.
#' @export
#'
#' @examples
#' \dontrun{
#' cox_all_genes(expression_vector_disease,time_vector,event_vector)
#' }
cox_all_genes <- function(eData,time_vector,event_vector){
  pb <- utils::txtProgressBar(min = 0, max = nrow(eData), style = 3)
  list_out <- list()
  for(i in 1:nrow(eData)){
    utils::setTxtProgressBar(pb, i)
    temp <- summary(survival::coxph(survival::Surv(time_vector,as.numeric(event_vector))~eData[i,]))$coefficients[1,]
    list_out[[i]] <- temp
  }
  df_out <- data.frame(do.call("rbind",list_out))
  colnames(df_out) <-  c("coef","exp_coef","se_coef","z","Pr_z")
  rownames(df_out) <-rownames(eData)
  df_out <- as.matrix(df_out)
  return(df_out)
}



#' Gene selector based on association to survival.
#'
#' @param cox_all A matrix output from the cox_all_genes function
#' @param percent A two element vector indicating the bottom and top quantiles for gene selection.
#'
#' @return Returns a list of genes that present z-values above or below the selected quantile thresholds.
#' @export
#'
#' @examples
#' \dontrun{
#' get_survival_related_genes(cox_all,percent)
#' }
get_survival_related_genes <- function(cox_all,percent = c(0.05,0.95)){
  genes_asso_surv <- rownames(cox_all[cox_all[,"z"] < stats::quantile(cox_all[,"z"],probs = percent[1]) | cox_all[,"z"] > stats::quantile(cox_all[,"z"],probs = percent[2]),])
  return(genes_asso_surv)

}
