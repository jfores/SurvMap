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
#' @param exp_to_res  An exponent in the form 1/n to wich the node sizes have to be exponentiatent in order to resize them.
#' @return Plots an interactive network using the visNetwork function.
#' @export
#'
#' @examples
#' \dontrun{
#' plot_mapper(mapper_list)}
plot_mapper <- function(mapper_list,trans_node_size = TRUE,exp_to_res = 1/2){
  arr_ind <- base::which(arr.ind = TRUE,mapper_list$adj_matrix == 1)
  df_out <- base::data.frame(base::rownames(mapper_list$adj_matrix)[arr_ind[,1]],base::colnames(mapper_list$adj_matrix)[arr_ind[,2]])
  df_out <- base::cbind(arr_ind,df_out)
  base::rownames(df_out) <- 1:base::nrow(df_out)
  base::colnames(df_out) <- c("from","to","from_Name","to_Name")
  nodes_to_net <- base::unique(base::data.frame(c(df_out[,1]-1,df_out[,2]-1),c(df_out[,3],df_out[,4])))
  nodes_to_net$node_size <- mapper_list$node_sizes
  if(trans_node_size){
    nodes_to_net$node_size <- (nodes_to_net$node_size)^exp_to_res
  }
  base::colnames(nodes_to_net) <- c("id","label","size")
  nodes_to_net$color <- map_to_color(base::log2(base::unlist(mapper_list$node_av_filt) + 2))
  edges_to_net <- df_out[,c(1,2)]-1
  base::colnames(edges_to_net) <- c("from","to")
  visNetwork::visNetwork(nodes_to_net,edges_to_net[!edges_to_net$from == edges_to_net$to,],)
}


#' Function to get colors like in ggplot2
#'
#' @param n number of colors.
#' @param h range.
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' ggplotColours(n,h)}
ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}


#' plot_genes_by_groups
#'
#' Generate plots of genes by groups.
#'
#' @param exp_data expression data.
#' @param genes genes to be plotted
#' @param out_perform_wil differential expression analysis results for multiple groups.
#' @param out_dir output directory.
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' plot_genes_by_groups(exp_data,genes,out_perform_wil)}
plot_genes_by_groups <- function(exp_data,genes,out_perform_wil){
  samp_in_nodes <- out_perform_wil[[2]]
  list_out <- list()
  for(i in 1:length(samp_in_nodes)){
    a_set <- rep(names(samp_in_nodes)[i],length(samp_in_nodes[[i]]))
    b_set <- samp_in_nodes[[i]]
    df_temp <- data.frame(a_set,b_set)
    list_out[[i]] <- df_temp
  }
  df_out <- do.call("rbind",list_out)
  print(df_out)
  group_data_ord <- unique(df_out[,1])[order(unique(df_out[,1]))]
  group_colors <- ggplotColours(length(group_data_ord))
  names(group_colors) <- group_data_ord
  #print(group_colors)
  exp_filt <- exp_data[,df_out[,2]]
  list_plots <- list()
  for(i in 1:length(genes)){
    df_temp <- data.frame(exp_filt[genes[i],],df_out[,1])
    colnames(df_temp) <- c("expression","group")
    p <- ggplot2::ggplot(df_temp, ggplot2::aes(x=group, y=expression,color = group)) +
      ggplot2::geom_violin() +  ggplot2::geom_boxplot(width=0.1, fill="white") + ggplot2::theme_classic() + ggplot2::ggtitle(paste("Gene: ",genes[i])) + ggplot2::scale_fill_manual(values=group_colors)
    list_plots[[i]] <- p
  }
  return(list_plots)
}


#' plot_heatmap_data
#'
#' @param out_wilc object returned by the differential gene expresson function.
#' @param selected_genes top differentially expressed genes.
#' @param exp_data gene expression matrix.
#' @param row_text_size size of the row labels.
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' plot_heatmap_data(out_wilc,selected_genes,exp_data,row_text_size = 10)}
plot_heatmap_data <- function(out_wilc,selected_genes,exp_data,row_text_size = 10){

  #Create samples in nodes data.
  list_out <- list()
  samp_in_nodes <- out_wilc[[2]]
  for(i in 1:length(samp_in_nodes)){
    a_set <- rep(names(samp_in_nodes)[i],length(samp_in_nodes[[i]]))
    b_set <- samp_in_nodes[[i]]
    df_temp <- data.frame(a_set,b_set)
    list_out[[i]] <- df_temp
  }
  df_out <- do.call("rbind",list_out)


  # Filter expression matrix.

  exp_data_filt <- exp_data[,df_out[,2]]

  # Getting group colors.

  group_data_ord <- unique(df_out[,1])[order(unique(df_out[,1]))]
  group_colors <- ggplotColours(length(group_data_ord))
  names(group_colors) <- group_data_ord

  # Create colour dataframe

  df_colors <- data.frame(names(group_colors),group_colors)
  colnames(df_colors) <- c("Nodes","colors")
  df_merged <- merge(df_out,df_colors,by.x = 1,by.y = 1,all.x = TRUE)
  print(head(df_merged))
  ha_data <- df_merged[,3]
  names(ha_data) <- df_merged[,1]

  # Filter and rename expression matgrix.

  exp_data_filt <- exp_data_filt[,df_merged[,2]]
  colnames(exp_data_filt) <- df_merged[,1]

  # Create ha object.

  ha_data <- group_colors
  ha = ComplexHeatmap::HeatmapAnnotation(bar = df_merged[,1], col = list(bar = group_colors))

  # Create color ramp.

  col_fun = circlize::colorRamp2(c(-4, 0, 4), c("red", "black", "green"))

  # Draw the heatmap.

  ComplexHeatmap::draw(ComplexHeatmap::Heatmap(t(scale(t(scale(exp_data_filt[selected_genes,])))),cluster_columns = F,col = col_fun,top_annotation = ha,cluster_rows = F,row_names_gp = gpar(fontsize = row_text_size),column_names_gp = gpar(fontsize = 0)))
}




