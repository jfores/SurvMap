#' map_to_color
#'
#' Auxiliary function that maps a numeric vector, the average node filtering function values, to a color vector.
#'
#' @param x A vector of numeric values storing the average filtering function values found in the samples placed into a specific node.
#' @param limits A two element numeric vector including the range of values. This is optional.
#'
#' @return A vector of the same length of x with colors ranging from blue to red.
#' @export
#'
#' @examples
#' \dontrun{
#' map_to_color(mapper_list)}
map_to_color <- function(x,limits=NULL){
  pallette_ob <-  grDevices::colorRampPalette(colors = c("blue","red"))(100)
  if(is.null(limits)){
    limits=range(x)}
  map_to_col <- pallette_ob[base::findInterval(x,base::seq(limits[1],limits[2],length.out=length(pallette_ob)+1), all.inside=TRUE)]
  return(map_to_col)
  }


#' plot_mapper
#'
#' This function produces an interactive network plot using the visNetork function.
#'
#' @param mapper_list A list produced as an output of the one_D_Mapper function.
#'
#' @return Plots an interactive network using the visNetwork function.
#' @export
#'
#' @examples
#' \dontrun{
#' plot_mapper(mapper_list)}
plot_mapper <- function(mapper_list){
  arr_ind <- base::which(arr.ind = TRUE,mapper_list$adj_matrix == 1)
  df_out <- base::data.frame(base::rownames(mapper_list$adj_matrix)[arr_ind[,1]],base::colnames(mapper_list$adj_matrix)[arr_ind[,2]])
  df_out <- base::cbind(arr_ind,df_out)
  base::rownames(df_out) <- 1:base::nrow(df_out)
  base::colnames(df_out) <- c("from","to","from_Name","to_Name")
  nodes_to_net <- base::unique(base::data.frame(c(df_out[,1]-1,df_out[,2]-1),c(df_out[,3],df_out[,4])))
  nodes_to_net$node_size <- mapper_list$node_sizes
  base::colnames(nodes_to_net) <- c("id","label","size")
  nodes_to_net$color <- map_to_color(base::log2(base::unlist(mapper_list$node_av_filt) + 2))
  edges_to_net <- df_out[,c(1,2)]-1
  base::colnames(edges_to_net) <- c("from","to")
  visNetwork::visNetwork(nodes_to_net,edges_to_net[!edges_to_net$from == edges_to_net$to,],)
}
