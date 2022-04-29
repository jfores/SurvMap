#' compute_power_tables
#'
#' Compute power tables for each dataset.
#'
#' @param list_datasets spitted list of expression matrices.
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' compute_power_tables(dis_st_mod,filt_vector,n_int,p)
#' }
compute_power_tables <- function(list_datasets){
  powers = c(seq(2, 20, by=1))
  powerTables = vector(mode = "list", length = length(list_datasets))
  for (set in 1:length(list_datasets)){
    powerTables[[set]] = list(data = WGCNA::pickSoftThreshold(list_datasets[[set]], powerVector=powers, corFnc = "bicor", networkType = "signed hybrid", blockSize = 10000, verbose = 2 )[[2]])
  }
  names(powerTables) <- names(list_datasets)
  return(powerTables)
}

#' plot_scale_free_fit
#'
#' Plots the scale-free topology fit for the different tested powers.
#'
#' @param pwer_table a single table of powers computed using the compute_power_tables function.
#' @param name_ds name of the dataset from which the power table derives.
#' @param cex_val cex value to control the size of the numbers printed in the plot.
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' plot_scale_free_fit(pwer_table,name_ds,cex_val = 0.7)
#' }
plot_scale_free_fit <- function(pwer_table,name_ds,cex_val = 0.7){
  plot(pwer_table$data[,1],-sign(pwer_table$data[,3])*pwer_table$data[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type = "n",ylim = c(0,1),main = paste("Scale independence in ", name_ds),cex.main = 0.7,cex.lab = 0.7)
  text(pwer_table$data[,1],-sign(pwer_table$data[,3])*pwer_table$data[,2],labels=pwer_table$data[,1],col="red",cex = cex_val)
  abline(h=0.90,col="red")
}

#' plot_mean_connectivity
#'
#' Plots the mean connectivity values for the different tested powers.
#'
#' @param pwer_table a single table of powers computed using the compute_power_tables function.
#' @param name_ds name of the dataset from which the power table derives.
#' @param cex_val cex value to control the size of the numbers printed in the plot.
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' plot_mean_connectivity(pwer_table,name_ds,cex_val = 0.7)
#' }
plot_mean_connectivity <- function(pwer_table,name_ds,cex_val = 0.7){
  plot(pwer_table$data[,1],pwer_table$data[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity",type = "n",ylim = c(0,max(pwer_table$data[,5])),main = paste("Mean connectivity in", name_ds),cex.main = 0.7,cex.lab = 0.7)
  text(pwer_table$data[,1],pwer_table$data[,5],labels=pwer_table$data[,1],col="red")
  abline(h=0.90,col="red")
}

#' plot_multiple_graphs
#'
#' plots scale-free topology fit or mean connectivity against power for multiple datasets using the output from the compute_power_tables function.
#'
#' @param power_tables output from the compute_power_tables function.
#' @param type fit or connectivity
#' @param file_to_save path and file to save the pdf.
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' plot_multiple_graphs(power_tables,type = c("fit","connectivity"),file_to_save)
#' }
plot_multiple_graphs <- function(power_tables,type = c("fit","connectivity"),file_to_save){
  sizeGrWindow(20,20)
  pdf(file =file_to_save,width = 12,height = 12)
  par(mfrow = c(ceiling(length(power_tables)/4),4));
  par(mar = c(3,3,1.5, 0.5))
  par(mgp = c(1.6, 0.6, 0))
  if(type == "fit"){
    for(i in 1:length(power_tables)){
      plot_scale_free_fit(power_tables[[i]],names(power_tables)[[i]])
    }
  }
  else if(type == "connectivity"){
    for(i in 1:length(power_tables)){
      plot_mean_connectivity(power_tables[[i]],names(power_tables)[[i]])
    }
  }
  dev.off()
}

#' select_powers_aux
#'
#' auxiliar function for the select powers function.
#'
#' @param pow_tab power table
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' select_powers_aux(pow_tab)
#' }
select_powers_aux <- function(pow_tab){
  selected_power <- pow_tab$data[-sign(pow_tab$data[,3])*pow_tab$data[,2] > 0.9,][1,1]
  if(is.na(selected_power)){
    selected_power <-  pow_tab$data[which.max(-sign(pow_tab$data[,3])*pow_tab$data[,2]),1]
  }
  return(selected_power)
}

#' select_powers
#'
#' Select powers that produce the best scale-free topology fit.
#'
#' @param power_tables power tables list produced by the compute_power_tables function.
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' select_powers(power_tables)
#' }
select_powers <- function(power_tables){
  powers <- unlist(lapply(power_tables,select_powers_aux))
  return(powers)
}

#' block_wise_mod
#'
#' Computes blockwiseConsensusModules
#'
#' @param list_exp list containing gene expressoin datasets.
#' @param powers power derived from the scale-free topology fit analysis.
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' block_wise_mod(list_exp, powers)
#' }
block_wise_mod <- function(list_exp, powers, mergeCutHeight = 0.25,deepSplit = 3,consensusQuantile = 0.25){
  WGCNA::allowWGCNAThreads()
  multiExpr = list()
  for( i in 1:length(list_exp)){
    multiExpr[[i]]<- list(data = list_exp[[i]])
  }
  mods <- WGCNA::blockwiseConsensusModules(multiExpr, maxBlockSize = 30000, corType = "bicor", power = powers, networkType = "signed hybrid", TOMDenom = "mean", checkMissingData = FALSE, deepSplit = deepSplit, reassignThresholdPS = 0, pamRespectsDendro = FALSE, mergeCutHeight = mergeCutHeight, numericLabels = TRUE, getTOMScalingSamples = TRUE, consensusQuantile = consensusQuantile, verbose = 3, indent = 2)
  return(mods)
}

#' meta_norm_nobs
#'
#' @param list_z_scores list of z_scores.
#' @param list_nobs list of number of observations.
#' @param column column
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' meta_norm_nobs(list_z_scores,list_nobs,column)
#' }
meta_norm_nobs <- function(list_z_scores,list_nobs,column){
  reduced_data <- Reduce("cbind",lapply(list_z_scores, function(x) x[,column]))
  reduced_data_nobs <- Reduce("cbind",lapply(list_nobs, function(x) x[,column]))
  bool_filt <- !(colSums(is.na(reduced_data)) > 0)
  reduced_data <- reduced_data[,bool_filt]
  reduced_data_nobs <- reduced_data_nobs[,bool_filt]
  reduced_data_nobs <- unname(reduced_data_nobs[1,])
  weight_st <- sqrt(reduced_data_nobs)
  number_of_studies <- ncol(reduced_data)
  z_comb_num <- rowSums(t(t(reduced_data) * weight_st))
  z_comb_den <- sqrt(sum(weight_st^2 ))
  z_comb <- z_comb_num / z_comb_den
  p_comb <- 2*pnorm(abs(z_comb), lower.tail = FALSE)
  return(list(z_comb,p_comb))
}

#' prepare_data_stouffer_p_val
#'
#' @param p_Data_split p_Data_split
#' @param ME_sig_Z ME_sig_Z
#' @param ME_sig_Nobs ME_sig_Nobs
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' prepare_data_stouffer_p_val(p_Data_split,ME_sig_Z,ME_sig_Nobs)
#' }
prepare_data_stouffer_p_val <- function(p_Data_split,ME_sig_Z,ME_sig_Nobs){
  list_z <- list()
  list_pval <- list()
  for(i in 1:(ncol(p_Data_split[[1]]))){
    print(i)
    temp_data <- meta_norm_nobs(ME_sig_Z,ME_sig_Nobs,i)
    list_z[[i]] <- temp_data[[1]]
    list_pval[[i]] <- temp_data[[2]]
  }
  list_z_df <- data.frame(do.call("cbind",list_z))
  colnames(list_z_df) <- colnames(p_Data_split[[1]])
  list_pval_df <- data.frame(do.call("cbind",list_pval))
  colnames(list_pval_df) <- colnames(p_Data_split[[1]])
  return(list(list_z_df,list_pval_df))
}

#' create_data_frame_from_multi
#'
#' @param vec_to_trans vector of covariate with more than 2 values.
#' @param row_names rownames of the original dataframe
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' create_data_frame_from_multi(vec_to_trans,row_names)
#' }
create_data_frame_from_multi <- function(vec_to_trans,row_names){
  values_vec <- unique(vec_to_trans)[order(unique(vec_to_trans))]
  df_data <- data.frame(matrix(0,nrow = length(vec_to_trans),ncol = length(values_vec)))
  for(i in 1:length(values_vec)){
    df_data[vec_to_trans %in% values_vec[i],i] <- 1
  }
  colnames(df_data) <- values_vec
  rownames(df_data) <- row_names
  return(df_data)
}

#' create_multi_expressoin_object
#'
#' @param test_splitted test_splitted
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' create_multi_expressoin_object(test_splitted)
#' }
create_multi_expressoin_object <- function(test_splitted){
  multiExpr = list()
  for( i in 1:length(test_splitted)){
    multiExpr[[i]]<- list(data = test_splitted[[i]])
  }
  return(multiExpr)
}

#' module_enrichment_analysis
#'
#' Performs enrichment analysis for WGCNA results.
#'
#' @param labels a vector containing the node to which each gene is assigend.
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' module_enrichment_analysis(tlabels)
#' }
module_enrichment_analysis <- function(labels,type_of_test = c("weight01","classic")){
  require(topGO)
  list_enrichment_out <- list()
  # Get background genes.
  bg_genes <- names(labels)
  # Get unique modules.
  vec_mods <- unique(labels)[order(unique(labels))]
  # Retrieve GO annotations from biomart
  db= biomaRt::useMart('ENSEMBL_MART_ENSEMBL',dataset='hsapiens_gene_ensembl',host="https://www.ensembl.org")
  go_ids= biomaRt::getBM(attributes=c('go_id', 'external_gene_name', 'namespace_1003'), filters='external_gene_name',values=bg_genes, mart=db)
  # Carry out analysis for each node.
  gene_2_GO=unstack(go_ids[,c(1,2)])
  for(i in 1:length(vec_mods)){
    print(paste("Carrying out enrichment for module: ",i, sep=""))
    candidate_list <- names(labels[labels == vec_mods[i]])
    keep = candidate_list %in% go_ids[,2]
    keep =which(keep==TRUE)
    candidate_list=candidate_list[keep]
    geneList=factor(as.integer(bg_genes %in% candidate_list))
    names(geneList)= bg_genes
    # Generating topGO object
    GOdata= new('topGOdata', ontology='BP', allGenes = geneList, annot = topGO::annFUN.gene2GO, gene2GO = gene_2_GO,nodeSize = 10)
    # Running topGO
    weight_fisher_result=topGO::runTest(GOdata, algorithm=type_of_test, statistic='fisher')
    # Formatting results table.
    allGO=usedGO(GOdata)
    all_res=topGO::GenTable(GOdata, weightFisher=weight_fisher_result, orderBy='weightFisher', topNodes=length(allGO))
    p.adj=p.adjust(all_res$weightFisher,method="BH")
    # create the file with all the statistics from GO analysis
    all_res_final=cbind(all_res,p.adj)
    all_res_final=all_res_final[order(all_res_final$p.adj),]
    list_enrichment_out[[i]] <- all_res_final
  }
  names(list_enrichment_out) <- paste("ME",vec_mods,sep = "")
  return(list_enrichment_out)
}

#' filt_fun_enrich
#'
#' Filter module enrichment analysis results.
#'
#' @param x enrichment list output.
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' filt_fun_enrich(x)
#' }
filt_fun_enrich <- function(x){
  selected <- c()
  for(i in 1:length(x)){
    print(i)
    try({
      if(class(x[[i]]) == "data.frame" & min(x[[i]][,"p.adj"],na.rm = TRUE) < 0.05){
        selected <- c(selected,i)
      }
    })
  }
  return(x[selected])
}
