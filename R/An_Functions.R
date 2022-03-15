#' get_information_from_results
#'
#' This functions retrieves information about the mapper results.
#'
#' @param res_mapper Result object from a mapper function.
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' get_information_from_results(res_mapper)
#' }
get_information_from_results <- function(res_mapper){
  n_nodes <- length(res_mapper$node_sizes)
  av_node_size <- mean(res_mapper$node_sizes)
  sd_node_size <- sd(res_mapper$node_sizes)
  n_connections <- sum(res_mapper$adj_matrix[upper.tri(res_mapper$adj_matrix)] == 1)
  prop_connections <- n_connections/length(res_mapper$adj_matrix[upper.tri(res_mapper$adj_matrix)])
  adj_mat <- res_mapper$adj_matrix
  adj_mat[lower.tri(adj_mat)] <- t(adj_mat)[lower.tri(adj_mat)]
  diag(adj_mat) <- 0
  ramifications_n <- colSums(adj_mat)-2
  ramifications_n[ramifications_n %in% c(-1,-2)] <- 0
  ramifications_n <- sum(ramifications_n)
  list_out <- list(n_nodes,av_node_size,sd_node_size,n_connections,prop_connections,ramifications_n)
  return(list_out)
}

#' generate_not_ov_clust
#'
#' Auxiliary function for the computation of principal components.
#'
#' @param out_one_D Results object obtained from the
#' @param exp_filt filtered expression matrix.
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' generate_not_ov_clust(out_one_D,exp_filt)
#' }
generate_not_ov_clust <- function(out_one_D,exp_filt){
  samples_in_nodes <- out_one_D$node_samples
  list_pc1 <- list()
  for(i in 1:length(samples_in_nodes)){
    #print(i)
    temp_samps <- samples_in_nodes[[i]]
    if(length(temp_samps) == 1){
      pc_1 <- exp_filt[,temp_samps]
      list_pc1[[i]] <- pc_1
    }else if(length(temp_samps) > 1){
      pc_1 <- stats::prcomp(t(exp_filt[,temp_samps]))
      pc_1 <- pc_1$rotation[,1]
      list_pc1[[i]] <- pc_1
    }
  }
  list_pc1 <- do.call("cbind",list_pc1)
  colnames(list_pc1) <- names(samples_in_nodes)
  return(list_pc1)
}




#' generate_df_samp_to_node
#'
#' Auxiliary function for the generation of the output data.frame.
#'
#' @param out_one_D Results object obtained from the
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' generate_df_samp_to_node(out_one_D)
#' }
generate_df_samp_to_node <- function(out_one_D){
  samples_in_nodes <- out_one_D$node_samples
  struc_out <- list()
  for(i in 1:length(samples_in_nodes)){
    df_temp <- data.frame(rep(names(samples_in_nodes)[i],length(samples_in_nodes[[i]])),samples_in_nodes[[i]])
    colnames(df_temp) <- c("Node","Sample")
    struc_out[[i]] <- df_temp
  }
  struc_out <- do.call("rbind",struc_out)
  return(struc_out)
}


#' create_output_df_unique_node
#'
#' Assigning samples to individual clusters.
#'
#' @param out_one_D Results object obtained from the
#' @param exp_filt filtered expression matrix
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' create_output_df_unique_node(out_one_D,exp_filt)
#' }
create_output_df_unique_node <- function(out_one_D,exp_filt){
  principal_components <- generate_not_ov_clust(out_one_D,exp_filt)
  struct_out <- generate_df_samp_to_node(out_one_D)
  print(colnames(struct_out))
  final_out <- data.frame(unique(struct_out[,2]),rep(NA,length(unique(struct_out[,2]))))
  colnames(final_out) <- c("sample","unique_cluster")
  for(i in 1:nrow(final_out)){
    #print(i)
    sample_in_nodes <- struct_out[struct_out[,2] %in% final_out[i,1],,drop = FALSE]
    if(nrow(sample_in_nodes) == 1){
      final_out[i,2] <- sample_in_nodes[,1]
    }else{
      cor_vec <- c()
      for(j in 1:nrow(sample_in_nodes)){
        cor_vec <- c(cor_vec,cor(exp_filt[,final_out[i,1]],principal_components[,sample_in_nodes[j,1]]))
        #cor_vec <- data.frame(cor_vec)
      }
      df_cor <- data.frame(sample_in_nodes[,1],cor_vec)
      colnames(df_cor) <- c("Node","Cor")
      selected <- df_cor[which.max(df_cor[,2]),1]
      final_out[i,2] <- selected
      #print(selected)
    }
  }
  return(final_out)
}

