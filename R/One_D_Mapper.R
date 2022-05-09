#' Extract intervals from filter function output values.
#'
#'This function constructs the intervals for a given filtering function.
#'
#' @param dis_st_mod Disease-state model matrix generated using the complete dataset.
#' @param filt_vector Vector obtained after applying the filtering function to the disease-state model data.
#' @param n_int Number of intervals to divide the filtering function in.
#' @param p Percentage of overlap between intervals.
#'
#' @return Returns the set of intervals for the filtering function.
#' @export
#'
#' @examples
#' \dontrun{
#' get_intervals_One_D(dis_st_mod,filt_vector,n_int,p)}
get_intervals_One_D <- function(dis_st_mod,filt_vector,n_int,p){
  range_filt <- max(filt_vector) - min(filt_vector)
  n_ov <- n_int -1
  l_int <- range_filt/(n_int - (n_ov*p))
  p_int <- p*l_int
  list_int <- list()
  list_int[[1]] <- c(min(filt_vector)-0.1,min(filt_vector + l_int))
  names(list_int)[1] <- "Level_1"
  for(i in 2:n_int){
    if(i < n_int){
      list_int[[i]] <- c(list_int[[i-1]][2]-p_int,list_int[[i-1]][2]-p_int + l_int)
      names(list_int)[i] <- paste("Level_",i,sep ="")
    }
    else if(i == n_int){
      list_int[[i]] <- c(list_int[[i-1]][2]-p_int,list_int[[i-1]][2]-p_int + l_int + 0.1)
      names(list_int)[i] <- paste("Level_",i,sep ="")
    }
  }
  return(list_int)
}



#' clust_lev
#'
#' Get clusters for a particular data level.
#'
#' @param dis_est_mod_lev Disease estate data for the whole datasets and the selected genes.
#' @param distance_type Type of distance , correlation, euclidean...
#' @param optimal_clust_mode Method for selection optimal number of clusters. For the moment only the standard method has been implemented.
#' @param n_bins_clust Number of bins in order to select the optimal threshold for number of cluster computations.
#' @param level_name Name of the studied level.
#'
#' @return Returns the samples included in each cluster for the specific level analyzed.
#' @export
#'
#' @examples
#' \dontrun{
#' clust_lev(dis_est_mod_lev,distance_type,optimal_clust_mode,n_bins_clust,level_name)}
clust_lev <- function(dis_est_mod_lev,distance_type = "cor",clust_type = "hierarchical",linkage_type = "single",optimal_clust_mode = "standard",n_bins_clust = 10,level_name = "level_1"){
  distances <- c("cor","euclidean")
  test_distances <- pmatch(distance_type,distances)
  if(is.na(test_distances)){
    stop(paste("Invalid distance selected. Choose one of the folowing:", paste(distances,collapse = ", ")))
  }
  clust_types <- c("hierarchical","PAM")
  test_clust_types <- pmatch(clust_type,clust_types)
  if(is.na(test_clust_types)){
    stop(paste("Invalid clustering method selected. Choose one of the folowing:", paste(clust_types,collapse = ", ")))
  }
  linkage_types <- c("single","average","complete")
  test_linkage_types <- pmatch(linkage_type,linkage_types)
  if(is.na(test_linkage_types)){
    stop(paste("Invalid linkage method selected. Choose one of the folowing:", paste(linkage_types,collapse = ", ")))
  }
  optimal_clust_modes <- c("standard","silhouette")
  test_optimal_clust_modes <- pmatch(optimal_clust_mode,optimal_clust_modes)
  if(is.na(test_optimal_clust_modes)){
    stop(paste("Invalid optimal cluster number method selected. Choose one of the folowing:", paste(optimal_clust_modes,collapse = ", ")))
  }
  if(!(distance_type %in% c("cor","euclidean"))){
    print("Provide one of the specified distance types")
    return(NULL)
  }
  else if(distance_type == "cor"){
    level_dist <- stats::as.dist(1-stats::cor(dis_est_mod_lev))
  }else{
    level_dist <- stats::dist(base::t(dis_est_mod_lev),method = distance_type)
  }
  max_dist_lev <- base::max(level_dist)
  if(clust_type == "hierarchical"){
    print("It is hierarchical clustering...")
    level_hclust_out <- stats::hclust(level_dist,method = linkage_type)
  #}
  if(!(optimal_clust_mode %in% c("standard","silhouette"))){
    print("Provide one of the specified optimal clustering method types")
    return(NULL)
  }else if(optimal_clust_mode == "standard"){
    heights <- level_hclust_out$height
    breaks_for_bins <- base::seq(from=min(heights), to=max_dist_lev, by=(max_dist_lev - base::min(heights))/n_bins_clust)
    #print(breaks_for_bins)
    histogram <- graphics::hist(c(heights,max_dist_lev), breaks=breaks_for_bins, plot=FALSE)
    #plot(histogram)
    #plot(histogram)
    hist_gap <- (histogram$counts == 0)
    if(all(!hist_gap)){
      print("There is no gap... therefore only one cluster...")
      cluster_indices_level <- base::rep(1,ncol(dis_est_mod_lev))
      names(cluster_indices_level) <- base::colnames(dis_est_mod_lev)
      return(cluster_indices_level)
    }else{
      #print("There is a gap... therefore potentially multiple clusters...")
      #print(histogram$mids)
      #print(hist_gap)
      trheshold_value_a <- histogram$mids[min(which(hist_gap == TRUE))]
      trheshold_value_b <- histogram$mids[min(which(hist_gap == TRUE))-1]
      #trheshold_value <- (trheshold_value_a + trheshold_value_b)/2
      trheshold_value <- trheshold_value_a
      #print(trheshold_value)
      #print(paste("The threshold value is: ",base::round(trheshold_value,digits = 2),sep=""))
      #plot(level_hclust_out)
      cluster_indices_level <- base::as.vector(stats::cutree(level_hclust_out, h=trheshold_value))
      base::names(cluster_indices_level) <- base::colnames(dis_est_mod_lev)
      return(cluster_indices_level)
    }
  }else if(optimal_clust_mode == "silhouette"){
    max_dist_lev <- base::max(level_dist)
    level_hclust_out <- stats::hclust(level_dist,method=linkage_type)
    n_clust <- c()
    av_sil <- c()
    for(i in 2:(length(level_hclust_out$order)-1)){
      n_clust <- c(n_clust,i)
      test <- cluster::silhouette(stats::cutree(level_hclust_out,i),level_dist)
      av_sil <- c(av_sil,mean(test[,3]))
    }
    #print(av_sil)
    if(max(av_sil) >= 0.25){
      op_clust <- n_clust[which.max(av_sil)]
      cluster_indices_level <- stats::cutree(level_hclust_out,op_clust)
      return(cluster_indices_level)
    }else{
      cluster_indices_level <- stats::cutree(level_hclust_out,1)
      return(cluster_indices_level)
    }
  }
  }#Este es el que he puesto extra.
  else if(clust_type == "PAM"){
    print("PAM clustering")
    av_sil <- c()
    n_clust <- c()
    for(i in 1:(ncol(dis_est_mod_lev)-1)){
      temp_clust <- cluster::pam(x =level_dist,diss = TRUE,k = i)
      if(i == 1){
        av_sil <- c(av_sil,0)
        n_clust <- c(n_clust,1)
      }else{
        av_sil <-c(av_sil,mean(cluster::silhouette(temp_clust$clustering,level_dist)[,3]))
        n_clust <- c(n_clust,i)
      }
    }
    if(max(av_sil) >= 0.25){
      op_clust <- n_clust[which.max(av_sil)]
      cluster_indices_level <- cluster::pam(x =level_dist,diss = TRUE,k = op_clust)$clustering
      return(cluster_indices_level)
    }else{
      cluster_indices_level <- rep(1,ncol(dis_est_mod_lev))
      names(cluster_indices_level) <- colnames(dis_est_mod_lev)
      return(cluster_indices_level)
    }
  }
}


#' samples_in_levels
#'
#' This function returns a list of vectors containing the samples included at each level.
#'
#' @param int_data Filter function intervals. The output list produced by the get_intervals_One_D function.
#' @param filter_function A vector with the filtering function values for each included sample. It is the output of any of the available filtering functions.
#'
#' @return A list with the samples included in each of the levels.
#' @export
#'
#' @examples
#' \dontrun{
#' samples_in_levels(int_data,filter_function)}
samples_in_levels <- function(int_data,filter_function){
  return(lapply(int_data,function(x,y) names(which(y >= x[1] & y < x[2])),filter_function))
}


#' clust_all_levels
#'
#' @param Dis_Est_Mod Disease estate model data for the complete dataset.
#' @param samp_in_lev A list including the samples placed at each level. It is the output of the samples_in_levels function.
#' @param distance_type The distance type employed.
#' @param optimal_clust_mode The optimal cluster method selection.
#' @param n_bins_clust The number of bines for histogram threshold selection. Used in the standard method for optimal number of clusters selection.
#'
#' @return Returns a list including the clusters found at each level.
#' @export
#'
#' @examples
#' \dontrun{
#' samples_in_levels(Dis_Est_Mod,samp_in_lev)}
clust_all_levels <- function(Dis_Est_Mod,samp_in_lev,distance_type = "cor",clust_type = "hierarchical",linkage_type = "single",optimal_clust_mode = "standard",n_bins_clust = 10){
  distances <- c("cor","euclidean")
  test_distances <- pmatch(distance_type,distances)
  if(is.na(test_distances)){
    stop(paste("Invalid distance selected. Choose one of the folowing:", paste(distances,collapse = ", ")))
  }
  clust_types <- c("hierarchical","PAM")
  test_clust_types <- pmatch(clust_type,clust_types)
  if(is.na(test_clust_types)){
    stop(paste("Invalid clustering method selected. Choose one of the folowing:", paste(clust_types,collapse = ", ")))
  }
  linkage_types <- c("single","average","complete")
  test_linkage_types <- pmatch(linkage_type,linkage_types)
  if(is.na(test_linkage_types)){
    stop(paste("Invalid linkage method selected. Choose one of the folowing:", paste(linkage_types,collapse = ", ")))
  }
  optimal_clust_modes <- c("standard","silhouette")
  test_optimal_clust_modes <- pmatch(optimal_clust_mode,optimal_clust_modes)
  if(is.na(test_optimal_clust_modes)){
    stop(paste("Invalid optimal cluster number method selected. Choose one of the folowing:", paste(optimal_clust_modes,collapse = ", ")))
  }
  list_out <- base::list()
  for(i in 1:base::length(samp_in_lev)){
    if(length(samp_in_lev[[i]]) > 2){
      clust_level_temp <- clust_lev(Dis_Est_Mod[,samp_in_lev[[i]]],distance_type = distance_type,clust_type = clust_type,linkage_type = linkage_type, optimal_clust_mode = optimal_clust_mode, n_bins_clust = n_bins_clust, level_name = base::paste("Level",i,sep="_"))
    }else{
      if(base::length(samp_in_lev[[i]]) < 3 & base::length(samp_in_lev[[i]]) > 0){
        clust_level_temp <- base::rep(1,length(samp_in_lev[[i]]))
        base::names(clust_level_temp) <- samp_in_lev[[i]]
      }else if(length(samp_in_lev[[i]]) == 0){
        clust_level_temp <- NA
      }
    }
    list_out[[i]] <- clust_level_temp
  }
  base::names(list_out) <- base::names(samp_in_lev)
  return(list_out)
}




#' levels_to_nodes
#'
#' Extract the nodes information based on the level clustering data.
#'
#' @param clust_all_levels_list A list obtained from the clust all levels function.
#'
#' @return A list including the sample content of each detected node.
#' @export
#'
#' @examples
#' \dontrun{
#' levels_to_nodes(clust_all_levels_list)}
levels_to_nodes <- function(clust_all_levels_list){
  node_counter <- c()
  nodes_list <- list()
  node_counter <- 1
  for(i in 1:base::length(clust_all_levels_list)){
    if(!base::all(base::is.na(clust_all_levels_list[[i]]))){
      clusters <- base::unique(clust_all_levels_list[[i]])
      for(j in 1:base::length(clusters)){
        nodes_list[[base::paste("Node",node_counter,sep="_")]] <- base::names(clust_all_levels_list[[i]][clust_all_levels_list[[i]] == clusters[j]])
        node_counter <- node_counter + 1
      }
    }
  }
  return(nodes_list)
}



#' compute_node_adjacency
#'
#' Computes the adjacency matrix.
#'
#' @param nodes_list A list with the sample content of each node.
#'
#' @return Returns a matrix object that stores a 1 if there are shared samples in two given nodos and a 0 otherwise.
#' @export
#'
#' @examples
#' \dontrun{
#' compute_node_adjacency(nodes_list)}
compute_node_adjacency <- function(nodes_list){
  adj_matrix <- base::matrix(0,nrow = base::length(nodes_list),ncol = base::length(nodes_list))
  for(i in 1:(base::length(nodes_list))){
    for(j in i:(base::length(nodes_list))){
      if(length(base::intersect(nodes_list[[i]],nodes_list[[j]])) > 0){
        adj_matrix[i,j] <- 1
      }
    }
  }
  base::colnames(adj_matrix) <- base::names(nodes_list)
  base::rownames(adj_matrix) <- base::names(nodes_list)
  return(adj_matrix)
}


#' one_D_Mapper
#'
#' Wrapping function to carry out the complete process.
#'
#' @param Ds_for_an Disease-state model data for the complete dataset.
#' @param filter_function Filtering function output. It is a vector containing the filtering function values for each sample.
#' @param n_int Number of intervals used to create the first sample partition based on filtering values.
#' @param p Maximum overlap between intervals. Expressed as a fraction from zero to one.
#' @param distance_type Type of distance for the clustering procedure. Can be cor or euclidean.
#' @param optimal_clust_mode Method to find the optimal number of clusters. Can be the standard method or the silhouette-based method.
#' @param n_bins_clust Number of bing to generate the histogram employed by the standard optimal number of cluster finder method.
#'
#' @return Returns a list with different elements.
#' @export
#'
#' @examples
#' \dontrun{
#' one_D_Mapper(Ds_for_an,filter_function,n_int = 10,p = 0.3,distance_type = "cor",optimal_clust_mode =  "standard",n_bins_clust = 10)}
one_D_Mapper <- function(Ds_for_an,filter_function,n_int = 10,p = 0.3,distance_type = "cor",clust_type = "hierarchical",linkage_type = "single", optimal_clust_mode =  "standard",n_bins_clust = 10){
  distances <- c("cor","euclidean")
  test_distances <- pmatch(distance_type,distances)
  if(is.na(test_distances)){
    stop(paste("Invalid distance selected. Choose one of the folowing:", paste(distances,collapse = ", ")))
  }
  clust_types <- c("hierarchical","PAM")
  test_clust_types <- pmatch(clust_type,clust_types)
  if(is.na(test_clust_types)){
    stop(paste("Invalid clustering method selected. Choose one of the folowing:", paste(clust_types,collapse = ", ")))
  }
  linkage_types <- c("single","average","complete")
  test_linkage_types <- pmatch(linkage_type,linkage_types)
  if(is.na(test_linkage_types)){
    stop(paste("Invalid linkage method selected. Choose one of the folowing:", paste(linkage_types,collapse = ", ")))
  }
  optimal_clust_modes <- c("standard","silhouette")
  test_optimal_clust_modes <- pmatch(optimal_clust_mode,optimal_clust_modes)
  if(is.na(test_optimal_clust_modes)){
    stop(paste("Invalid optimal cluster number method selected. Choose one of the folowing:", paste(optimal_clust_modes,collapse = ", ")))
  }
  #Getting intervals.
  int_data <- get_intervals_One_D(Ds_for_an,filter_function,n_int = n_int,p = p)
  #Getting samples on each interval.
  sam_in_lev <- samples_in_levels(int_data,filter_function)
  #Clustering all levels.
  test_clust_all_levels <- clust_all_levels(Dis_Est_Mod = Ds_for_an,samp_in_lev = sam_in_lev,distance_type = distance_type, clust_type = clust_type, optimal_clust_mode =  optimal_clust_mode,n_bins_clust = n_bins_clust)
  #Transforming levels into nodes.
  node_samples <- levels_to_nodes(test_clust_all_levels)
  #Computing adjacency matrix.
  adj_matrix_out <- compute_node_adjacency(node_samples)
  #Creating parameter vector
  par_vec <- c(n_int,p,distance_type,optimal_clust_mode,n_bins_clust)
  names(par_vec) <- c("n_int","p","distance_type","optimal_clust_mode","n_bins_clust")
  data_out <- list()
  if(clust_type == "PAM"){
    par_vec[4] <- "silhouette"
    par_vec[5] <- NA
  }
  if(optimal_clust_mode == "silhouette"){
    par_vec[5] <- NA
  }
  #Generating the output data list.
  data_out$int_data <- int_data
  data_out$samp_in_lev <- sam_in_lev
  data_out$clust_all_levels <- test_clust_all_levels
  data_out$node_samples <- node_samples
  data_out$node_sizes <- unlist(lapply(node_samples,length))
  data_out$node_av_filt <- lapply(node_samples,function(x,y) mean(y[x]),filter_function)
  data_out$adj_matrix <- adj_matrix_out
  data_out$parameters <- par_vec
  return(data_out)
}

