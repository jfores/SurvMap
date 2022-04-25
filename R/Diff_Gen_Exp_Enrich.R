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
#' \dontrun{
#' select_top_diff_genes_groups(res_diff_exp,n_genes)
#' }
select_top_diff_genes_groups <- function(res_diff_exp,n_genes){
  selected_genes <- unique(unlist(lapply(res_diff_exp[[1]],aux_function_select_top_diff,n_genes)))
  return(selected_genes)
}


#' prepare_diff_meth_res
#'
#' Prepares differential methylation analysis results, and retrieves the top n probes presenting highest up- and -down differential methylation
#'
#' @param res_diff output object from perform_wilcoxon_each function.
#' @param gen_inf Information about the genes associated to each methylation probe obtained after applying the rewRanges function to a summarized experiment object
#' @param n_genes Number of top methylated genes to retrieve for each Node.
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' prepare_diff_meth_res(res_diff,gen_inf,n_genes)
#' }
prepare_diff_meth_res <- function(res_diff,gen_inf,n_genes){
  res_diff_one <- res_diff[[1]]
  res_diff_one_merged <- lapply(res_diff_one,merge,gen_inf,by.x = "row.names",by.y = "probeID" )
  out_res_top <- lapply(res_diff_one_merged,fun_to_app_prep,n_genes)
  return(out_res_top)
}


#' fun_to_app_prep
#'
#' Auxiliar function for prepare_diff_meth_res
#'
#' @param x A dataframe
#' @param n_genes Number of genes showing the highest differences in methylation.
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' fun_to_app_prep(x,n_genes)
#' }
fun_to_app_prep <- function(x,n_genes){

  x_up <- x[order(x[,4],decreasing = F),]
  x_up <- x_up[!x_up[,11] == "",]
  row_temp <- 1
  list_genes_temp <- c()
  list_selected_rows <- c()
  keep_searching <- TRUE
  while(keep_searching){
    print(row_temp)
    if(!x_up[row_temp,11] %in% list_genes_temp){
      list_selected_rows <- c(list_selected_rows,row_temp)
      list_genes_temp <- c(list_genes_temp,x_up[row_temp,11])
    }
    row_temp <- row_temp + 1
    if(length(list_genes_temp) >= n_genes){
      keep_searching <- FALSE
    }
  }
  print(list_genes_temp)
  x_up <- x_up[list_selected_rows,]

  x_down <- x[order(x[,5],decreasing = F),]
  x_down <- x_down[!x_down[,11] == "",]
  row_temp <- 1
  list_genes_temp <- c()
  list_selected_rows <- c()
  keep_searching <- TRUE
  while(keep_searching){
    print(row_temp)
    if(!x_down[row_temp,11] %in% list_genes_temp){
      list_selected_rows <- c(list_selected_rows,row_temp)
      list_genes_temp <- c(list_genes_temp,x_down[row_temp,11])
    }
    row_temp <- row_temp + 1
    if(length(list_genes_temp) >= n_genes){
      keep_searching <- FALSE
    }
  }
  x_down <- x_down[list_selected_rows,]
  print(list_genes_temp)
  return(rbind(x_up,x_down))
}


#' perform_GSVA_dataset
#'
#' Carries out GSVA for the complete dataset and camputes the averange values by pathways and group of samples.
#'
#' @param expression Gene expression matrix.
#' @param surv_Map_Out Output from SurvMap analyses.
#' @param path_to_gmt_file Pathway to the gtm file containing the pathways atabase.
#' @param thr_groups threshold of the number of samples used to filter out small nodes.
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' perform_GSVA_dataset(expression,surv_Map_Out,path_to_gmt_file,thr_groups = 20)
#' }
perform_GSVA_dataset <- function(expression,surv_Map_Out,path_to_gmt_file,thr_groups = 20){

  gene_sets_data <- GSA::GSA.read.gmt(path_to_gmt_file)
  gene_sets <- gene_sets_data$genesets
  names(gene_sets) <- gene_sets_data$geneset.names

  uniq_samp <- surv_Map_Out$Unique_Samp_Node
  rownames(uniq_samp) <- uniq_samp[,1]

  gsva_res <- GSVA::gsva(expression, gene_sets, min.sz=5, max.sz=500)
  print(dim(gsva_res))
  intersected_samples <- intersect(uniq_samp[,1],colnames(gsva_res))
  gsva_res <- gsva_res[,intersected_samples]
  uniq_samp <- uniq_samp[intersected_samples,]

  selected_nodes <- names(table(uniq_samp[,2]) > thr_groups)[table(uniq_samp[,2]) > thr_groups]
  uniq_samp <- uniq_samp[uniq_samp[,2] %in% selected_nodes,]
  print(table(uniq_samp[,2]))

  gsva_res <- gsva_res[,rownames(uniq_samp)]

  averaged_data <- t(apply(gsva_res,1,axuiliar_function_to_mean,y = uniq_samp[,2]))
  list_out <- list(gsva_res,averaged_data)
  return(list_out)
}

#' axuiliar_function_to_mean
#'
#' Auxiliar function that computes the mean of a vector based on a group.
#'
#' @param x vector of numeric values.
#' @param y vector indicating the group to which each numeric value belong.
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' axuiliar_function_to_mean(x,y)
#' }
axuiliar_function_to_mean <- function(x,y){
  return(tapply(x,y,mean))
}
