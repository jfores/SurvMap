


get_intervals_One_D <- function(dis_st_mod,filt_vector,n_int,p){
  range_filt <- max(filt_vector) - min(filt_vector)
  n_ov <- n_int -1
  l_int <- range_filt/(n_int - (n_ov*p))
  p_int <- p*l_int
  list_int <- list()
  list_int[[1]] <- c(min(filt_vector),min(filt_vector + l_int))
  for(i in 2:n_int){
    list_int[[i]] <- c(list_int[[i-1]][2]-p_int,list_int[[i-1]][2]-p_int + l_int)
  }
  return(list_int)
}
