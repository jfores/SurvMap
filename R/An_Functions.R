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
      pc_1 <- stats::prcomp(t(exp_filt[,temp_samps]),center = FALSE,scale. = TRUE)
      pc_1 <- pc_1$rotation[,1]
      av_exp_vec <- scale(apply(exp_filt[,temp_samps],1,mean),center = FALSE,scale = TRUE)
      if(cor(pc_1,av_exp_vec) < 0){
        pc_1 <- -pc_1
      }
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
        cor_vec <- c(cor_vec,cor(scale(exp_filt[,final_out[i,1]], center = FALSE, scale = TRUE),principal_components[,sample_in_nodes[j,1]]))
        #cor_vec <- data.frame(cor_vec)
      }
      df_cor <- data.frame(sample_in_nodes[,1],cor_vec)
      colnames(df_cor) <- c("Node","Cor")
      print(df_cor)
      selected <- df_cor[which.max(df_cor[,2]),1]
      final_out[i,2] <- selected
      #print(selected)
    }
  }
  return(final_out)
}



#' percen_table
#'
#' Computes percentages of categorical variables in groups of samples.
#'
#' @param vec_1 vector of groups
#' @param vec_2 vector of categorical variables.
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' create_output_df_unique_node(out_one_D,exp_filt)
#' }
percent_table <- function(vec_1,vec_2){
  table_vals <- table(vec_1,vec_2)
  print(table_vals)
  col_summed <- colSums(table_vals)
  print(col_summed)
  return(t(t(table_vals)/col_summed*100))
}

#' compute_u_mat_and_eigen
#'
#' Computes the matrix U and the singular values for a subset of samples placed in a node.
#'
#' @param dis_comp_node Subset of the disease component matrix.
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' compute_u_mat_and_eigen(dis_comp_node)
#' }
compute_u_mat_and_eigen <- function(dis_comp_node){
  C <- dis_comp_node%*%t(dis_comp_node)
  C <- (1/sum(diag(C)))*C
  base <- svd(C)$u
  eigen_values <- svd(C)$d
  return(list(base,eigen_values))
}

#' assign_promiscuous_samples
#'
#' Assigns promiscuous samples to unique nodes.
#'
#' @param out_one_D Output from one_D_Mapper function.
#' @param Dis_Comp_Data Disease component matrix.
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' assign_promiscuous_samples(out_one_D,Dis_Comp_Data)
#' }
assign_promiscuous_samples <- function(out_one_D,Dis_Comp_Data){
  samples_in_nodes <- out_one_D$node_samples
  list_Dis_Comp_Data_Node <- list()
  for(i in 1:length(samples_in_nodes)){
    list_Dis_Comp_Data_Node[[i]] <- Dis_Comp_Data[,samples_in_nodes[[i]]]
  }
  names(list_Dis_Comp_Data_Node) <- names(samples_in_nodes)
  u_mat_and_eigen_vals <- lapply(list_Dis_Comp_Data_Node,compute_u_mat_and_eigen)
  names(u_mat_and_eigen_vals) <- names(samples_in_nodes)
  struct_out <- generate_df_samp_to_node(out_one_D)
  final_out <- data.frame(unique(struct_out[,2]),rep(NA,length(unique(struct_out[,2]))))
  colnames(final_out) <- c("sample","unique_cluster")
  for(i in 1:nrow(final_out)){
    sample_to_nodes <- struct_out[struct_out[,2] %in% final_out[i,1],,drop = FALSE]
    if(nrow(sample_to_nodes) == 1){
      final_out[i,2] <- sample_to_nodes[,1]
    }else{
      prob_vectors <- c()
      for(j in 1:nrow(sample_to_nodes)){
        vector_to_test <- Dis_Comp_Data[,final_out[i,1]]
        print(sample_to_nodes[j,1])
        u_mat_temp <- u_mat_and_eigen_vals[[sample_to_nodes[j,1]]][[1]]
        eigen_vec_temp <- u_mat_and_eigen_vals[[sample_to_nodes[j,1]]][[2]]
        vector_to_test <- vector_to_test/sqrt(sum(vector_to_test^2))
        n_vector_to_test <- t(u_mat_temp)%*%vector_to_test
        prob <- sum(n_vector_to_test^2*eigen_vec_temp)
        print(prob)
        prob_vectors <- c(prob_vectors,prob)
      }
      df_probs <- data.frame(sample_to_nodes[,1],prob_vectors)
      colnames(df_probs) <- c("Node","Cor")
      selected <- df_probs[which.max(df_probs[,2]),1]
      final_out[i,2] <- selected
    }
  }
  return(final_out)
}



