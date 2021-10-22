

flatten_healthy <- function(normal_tissue_data){
  df_out <- normal_tissue_data
  for(i in 1:ncol(normal_tissue_data)){
    df_out[,i] <- fitted(lm(normal_tissue_data[,i] ~ 0 + ., data = data.frame(normal_tissue_data)[,-i]))
  }
  return(df_out)
}
