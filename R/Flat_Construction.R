

#' Computes lambda
#'
#'
#' Computes the value of lambda as defined in: The Optimal Hard Thresholdfor Singular Values is \deqn{\sqrt(4/ 3)}
#'
#' @param bet numeric
#'
#' @return numeric
#' @export
#'
#' @examples
#' get_lambda(0.3)
get_lambda <- function(bet){
  w = (8*bet) / ((bet + 1) + sqrt(bet^2 + 14*bet + 1))
  lambda_star <- sqrt(2*(bet + 1) + w)
  return(lambda_star)
}


fun_to_int <- function(t,bet){
  b_p <- (1 + sqrt(bet))^2
  b_m <- (1 - sqrt(bet))^2
  numerator <- sqrt((b_p - t)*(t - b_m))
  denominator <- 2*pi*t*bet
  res <- numerator/denominator
  return(res)
}


get_mu_beta <- function(bet){
  thresh_diff <- 1e-10
  lbond <- (1 - sqrt(bet))^2
  hbond <- (1 + sqrt(bet))^2
  seqvals <- seq(lbond,hbond,length.out = 100)
  end_loop <- FALSE
  counter <- 1
  while(end_loop == FALSE){
    #print(counter)
    values_int <- c()
    for (i in 1:length(seqvals)){
      res <- integrate(fun_to_int,lbond,seqvals[i],bet)$value
      values_int <- c(values_int,res)
    }
    if(abs(max(values_int[values_int < 0.5]) - 0.5) < thresh_diff & abs(min(values_int[values_int > 0.5]) - 0.5) < thresh_diff){
      #print("Done")
      final_seqval <- (max(seqvals[values_int < 0.5]) + min(seqvals[values_int > 0.5]))/2
      end_loop <- TRUE
    }else{
      seqvals <- seq(seqvals[max(which(values_int < 0.5))],seqvals[min(which(values_int > 0.5))],length.out = 100)
    }
    counter <- counter + 1
  }
  return(final_seqval)
}


get_omega <- function(bet){
  lamb <- get_lambda(bet)
  mu_beta <- get_mu_beta(bet)
  omega <- lamb/sqrt(mu_beta)
  return(omega)
}
