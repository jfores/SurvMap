



map_to_color <- function(x,limits=NULL){
  pallette_ob <-  grDevices::colorRampPalette(colors = c("blue","red"))(100)
  if(is.null(limits)){
    limits=range(x)}
  map_to_col <- pallette_ob[base::findInterval(x,base::seq(limits[1],limits[2],length.out=length(pallette_ob)+1), all.inside=TRUE)]
  return(map_to_col)
  }


plot_mapper <- function(mapper_list){
  NULL
}
