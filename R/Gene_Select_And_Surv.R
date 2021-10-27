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


get_survival_related_genes <- function(cox_all,percent = c(0.05,0.95)){
  genes_asso_surv <- rownames(cox_all[cox_all[,"z"] < stats::quantile(cox_all[,"z"],probs = percent[1]) | cox_all[,"z"] > stats::quantile(cox_all[,"z"],probs = percent[2]),])
  return(genes_asso_surv)

}
