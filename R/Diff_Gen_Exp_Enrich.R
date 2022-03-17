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
