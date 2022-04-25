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
#' compute_power_tables(dis_st_mod,filt_vector,n_int,p)}
compute_power_tables <- function(list_datasets){
  powers = c(seq(2, 20, by=1))
  powerTables = vector(mode = "list", length = length(test_splitted))
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
#' plot_scale_free_fit(pwer_table,name_ds,cex_val = 0.7)}
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
#' plot_mean_connectivity(pwer_table,name_ds,cex_val = 0.7)}
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
#' plot_multiple_graphs(power_tables,type = c("fit","connectivity"),file_to_save)}
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
#' select_powers_aux(pow_tab)}
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
#' select_powers(power_tables)}
select_powers <- function(power_tables){
  powers <- unlist(lapply(power_tables,select_powers_aux))
  return(powers)
}

