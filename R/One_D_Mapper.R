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
  list_int[[1]] <- c(min(filt_vector),min(filt_vector + l_int))
  names(list_int)[1] <- "Level_1"
  for(i in 2:n_int){
    list_int[[i]] <- c(list_int[[i-1]][2]-p_int,list_int[[i-1]][2]-p_int + l_int)
    names(list_int)[i] <- paste("Level_",i,sep ="")
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
#' clust_lev(dis_est_mod_lev,distance_type,optimal_clust_mode,n_bins_clust,level_name)}silhouette
clust_lev <- function(dis_est_mod_lev,distance_type = c("cor","euclidean"),optimal_clust_mode = c("standard",""),n_bins_clust = 10,level_name = "level_1"){
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
  level_hclust_out <- stats::hclust(level_dist,method="single")
  if(!(optimal_clust_mode %in% c("standard","silhouette"))){
    print("Provide one of the specified optimal clustering method types")
    return(NULL)
  }else if(optimal_clust_mode == "standard"){
    heights <- level_hclust_out$height
    breaks_for_bins <- base::seq(from=min(heights), to=max_dist_lev, by=(max_dist_lev - base::min(heights))/n_bins_clust)
    #print(breaks_for_bins)
    histogram <- graphics::hist(c(heights,max_dist_lev), breaks=breaks_for_bins, plot=FALSE)
    plot(histogram)
    #plot(histogram)
    hist_gap <- (histogram$counts == 0)
    if(all(!hist_gap)){
      print("There is no gap... therefore only one cluster...")
      cluster_indices_level <- base::rep(1,ncol(dis_est_mod_lev))
      names(cluster_indices_level) <- base::colnames(dis_est_mod_lev)
      return(cluster_indices_level)
    }else{
      print("There is a gap... therefore potentially multiple clusters...")
      print(histogram$mids)
      print(hist_gap)
      trheshold_value_a <- histogram$mids[min(which(hist_gap == TRUE))]
      trheshold_value_b <- histogram$mids[min(which(hist_gap == TRUE))-1]
      #trheshold_value <- (trheshold_value_a + trheshold_value_b)/2
      trheshold_value <- trheshold_value_b
      print(trheshold_value)
      print(paste("The threshold value is: ",base::round(trheshold_value,digits = 2),sep=""))
      #plot(level_hclust_out)
      cluster_indices_level <- base::as.vector(stats::cutree(level_hclust_out, h=trheshold_value))
      base::names(cluster_indices_level) <- base::colnames(dis_est_mod_lev)
      return(cluster_indices_level)
    }
  }else if(optimal_clust_mode == "silhouette"){
    print("To do...")
    max_dist_lev <- base::max(level_dist)
    level_hclust_out <- stats::hclust(level_dist,method="single")
    n_clust <- c()
    av_sil <- c()
    for(i in 2:(length(level_hclust_out$order)-1)){
      n_clust <- c(n_clust,i)
      test <- cluster::silhouette(stats::cutree(level_hclust_out,i),level_dist)
      av_sil <- c(av_sil,mean(test[,3]))
    }
    print(max(av_sil))
    if(max(av_sil) >= 0.25){
      op_clust <- n_clust[which.max(av_sil)]
      cluster_indices_level <- stats::cutree(level_hclust_out,op_clust)
      return(cluster_indices_level)
    }else{
      cluster_indices_level <- stats::cutree(level_hclust_out,1)
      return(cluster_indices_level)
    }
  }
}

