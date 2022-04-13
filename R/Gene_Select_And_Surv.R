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



#' gene_selection_surv
#'
#' Selects markers for mapper based on the product of standard deviation plus one times the Z score obtained thorugh fitting a cox proportional hazard model to the expression of each gene.
#'
#' @param D_Comp Disease component object
#' @param p_Data Data.frame with the phenotype data.
#' @param Status_Col_Name Column for sample filtering.
#' @param Status_Value Value for the sample filtering co-variate.
#' @param Cox_All Output from the cox_all_genes function.
#' @param n_top Number of markers presenting the maximum values for the product.
#' @param type_sel Select top from bottom and top or from absolute value.
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' gene_selection_surv(D_Comp,p_Data,Status_Col_Name,Status_Value,Cox_All,n_top)
#' }
gene_selection_surv <- function(D_Comp,p_Data,Status_Col_Name,Status_Value,Cox_All,n_top,type_sel = c("Top_Bot","Abs")){
  if(type_sel == "Top_Bot"){
    print("Is Top_Bot")
    probes_test <- (apply(D_Comp[,p_Data[,Status_Col_Name] == Status_Value],1,sd)+1) * Cox_All[,4]
    if(n_top %% 2 == 0){
      n_top <- n_top / 2
      print(n_top)
    }else{
      n_top <- (n_top + 1)/2
      print(n_top)
    }
    selected_probes <- names(c(probes_test[order(probes_test,decreasing = T)][1:n_top],probes_test[order(probes_test,decreasing = F)][1:n_top]))
    return(selected_probes)
  }else if(type_sel == "Abs"){
    print("Is Abs")
    probes_test <- (apply(D_Comp[,p_Data[,Status_Col_Name] == Status_Value],1,sd)+1) * abs(Cox_All[,4])
    selected_probes <- names(probes_test[order(probes_test,decreasing = F)][1:n_top])
    return(selected_probes)
  }
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


#' surivival_analysis_multiple_groups
#'
#' Survival analysis and plottding of the clustering results. https://www.reneshbedre.com/blog/survival-analysis.html
#'
#' @param pheno_data data.frame with phenotype data. Must include a columns called pCh_DFS_E, pCh_DFS_T, and pCh_Status including the event, time to relapse, and the status of the sample T: Disease NT: healthy
#' @param out_one_D Mapper output.
#' @param thr_groups Minimum number of samples to include a node in the analysis.
#' @param ylim_val y-axis limits for survival plot.
#' @param xlim_val x-axis limits for survival plot.
#' @param type Divide by node or by pam50 in the latter case a column pam50_frma containing the pam50 classification will be requested.
#' @param selected_nodes_2 subset of nodes to compare.
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' surivival_analysis_multiple_groups(pheno_data,out_one_D,thr_groups = 50,ylim_val = c(0.5, 1),xlim_val =  c(0,200),type = "node")
#' }
surivival_analysis_multiple_groups <- function(pheno_data,out_one_D,thr_groups = 50,ylim_val = c(0.5, 1),xlim_val =  c(0,200),type = "node",selected_nodes_2 = ""){
  univoq_group <- out_one_D$Unique_Samp_Node
  p_merged <- merge(pheno_data,univoq_group,by.x = 1,by.y = 1 )
  p_merged <- p_merged[p_merged$pCh_Status == "T",]
  p_merged <- p_merged[!(is.na(p_merged$pCh_DFS_E) | is.na(p_merged$pCh_DFS_T)),]
  selected_nodes <- names(table(p_merged$unique_cluster) > thr_groups)[table(p_merged$unique_cluster) > thr_groups]
  p_merged <- p_merged[p_merged$unique_cluster %in% selected_nodes,]
  p_merged$pCh_DFS_T <- as.numeric(p_merged$pCh_DFS_T)
  p_merged$pCh_DFS_E <- as.numeric(p_merged$pCh_DFS_E)
  p_merged <- p_merged
  group_data_ord <- unique(p_merged$unique_cluster)[order(unique(p_merged$unique_cluster))]
  group_colors <- ggplotColours(length(group_data_ord))
  names(group_colors) <- group_data_ord
  print(group_colors)
  if(! c("") %in%selected_nodes_2){
    p_merged <- p_merged[p_merged$unique_cluster %in% selected_nodes_2,]
    group_colors <- group_colors[names(group_colors) %in% selected_nodes_2]
    print(group_colors)
  }
  p_merged <<- p_merged
  surv = survival::Surv(time = as.numeric(p_merged$pCh_DFS_T), event = as.numeric(p_merged$pCh_DFS_E))
  if(type == "node"){
    fit <- survival::survfit(survival::Surv(time = p_merged$pCh_DFS_T, event = p_merged$pCh_DFS_E)~p_merged$unique_cluster)
    log_rank_test <- survival::survdiff(formula = survival::Surv(time = as.numeric(p_merged$pCh_DFS_T), event = as.numeric(p_merged$pCh_DFS_E)) ~ p_merged$unique_cluster)
  }else if(type == "pam"){
    fit <- survival::survfit(survival::Surv(time = p_merged$pCh_DFS_T, event = p_merged$pCh_DFS_E)~p_merged$pam50_frma)
    log_rank_test <- survival::survdiff(formula = survival::Surv(time = as.numeric(p_merged$pCh_DFS_T), event = as.numeric(p_merged$pCh_DFS_E)) ~ p_merged$pam50_frma)
  }
  plot_out <- survminer::ggsurvplot(fit = fit,data = p_merged, pval = TRUE, surv.median.line = "hv", xlab = "Survival time", ylab = "Survival probability",ylim = ylim_val,xlim = xlim_val,palette = unname(group_colors))
  p_merged_temp <- p_merged
  rm(list = "p_merged",envir = globalenv())
  return(list(fit,p_merged_temp,surv,plot_out,log_rank_test))
}
