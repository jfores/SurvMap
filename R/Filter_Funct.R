

lp_norm_k_powers <- function(exp_matrix,p,k){
  if(is.matrix(exp_matrix)){
    return(apply(exp_matrix,2,function(x) (sum(abs(x)^p)^(k/p))))
  }else{
    print("A  matrix object should be provided.")
  }
}
