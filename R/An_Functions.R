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
