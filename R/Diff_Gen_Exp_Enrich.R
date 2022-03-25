#' perform_wilcoxon_each
#'
#' @param exp_data Gene expression data for the complete dataset.
#' @param pheno_data Pheno data with columns pCh_Sample_Name "sample names", pCh_Status "T for disease and NT for non-disease samples"
#' @param surv_map_res surv_map_res object
#' @param thr_groups Threshold to filter out groups or nodes with less than n samples.
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' perform_wilcoxon_each(exp_data,pheno_data,surv_map_res,thr_groups = 20)
#' }
perform_wilcoxon_each <- function(exp_data,pheno_data,surv_map_res,thr_groups = 20){
  univoq_group <- surv_map_res$Unique_Samp_Node
  p_merged <- merge(pheno_data,univoq_group,by.x = 1,by.y = 1 )
  p_merged <- p_merged[p_merged$pCh_Status == "T",]
  selected_nodes <- names(table(p_merged$unique_cluster) > thr_groups)[table(p_merged$unique_cluster) > thr_groups]
  p_merged <- p_merged[p_merged$unique_cluster %in% selected_nodes,]
  unique_groups <- unique(p_merged$unique_cluster)
  out_results <- list()
  list_samples_in_groups <- list()
  for(i in 1:length(unique_groups)){
    print(paste("Computing differential expression for group ",unique_groups[i],sep=" "))
    samples_in_group <- p_merged[p_merged$unique_cluster == unique_groups[i],"pCh_Sample_Name"]
    list_samples_in_groups[[i]] <- samples_in_group
    samples_not_in_group <- p_merged[!p_merged$pCh_Sample_Name %in% samples_in_group,"pCh_Sample_Name"]
    exp_data_filt <- exp_data[,c(samples_in_group,samples_not_in_group)]
    factor_to_comp <- c(rep(0,length(samples_in_group)),rep(1,length(samples_not_in_group)))
    out_results[[i]] <- GSALightning::wilcoxTest(exp_data_filt,as.factor(factor_to_comp))
  }
  names(out_results) <- unique_groups
  names(list_samples_in_groups) <- unique_groups
  return(list(out_results,list_samples_in_groups))
}


#' aux_function_select_top_diff
#'
#' Axuliar function to select_top_diff_genes_groups
#'
#' @param x res_diff_exp object
#' @param n_genes number of top genes.
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' aux_function_select_top_diff(x,n_genes)
#' }
aux_function_select_top_diff <- function(x,n_genes){
  return(c(rownames(x[order(x[,3]),])[1:n_genes],rownames(x[order(x[,4]),])[1:n_genes]))
}

#' select_top_diff_genes_groups
#'
#' Select top differentially expressed genes for each group.
#'
#' @param res_diff_exp res_diff_exp object.
#' @param n_genes top differentially expressed genes.
#'
#' @return
#' @export
#'
#' @examples
#' select_top_diff_genes_groups(res_diff_exp,n_genes)
#' }
select_top_diff_genes_groups <- function(res_diff_exp,n_genes){
  selected_genes <- unique(unlist(lapply(res_diff_exp[[1]],aux_function_select_top_diff,n_genes)))
  return(selected_genes)
}



