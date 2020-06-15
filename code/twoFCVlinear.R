twoFCVlinear <- function(n,Y,j0,qmf,D,wc,sigma_u){
  n1 <- n/2
  lD <- length(D)
  wc_odd_lin  <- rep(0,n1)
  wc_even_lin <- rep(0,n1)
  Y_odd <- Y[seq(1,n,2)]
  Y_even <- Y[seq(2,n,2)]

  wc_odd  <- FWT_PO(Y_odd,j0,qmf)
  wc_even <- FWT_PO(Y_even,j0,qmf)
  hat_s_odd  <- matrix(0,nrow=lD,ncol=n1)
  hat_s_even <- matrix(0,nrow=lD,ncol=n1)
  hat_s_m <- matrix(rep(0,lD*n),nrow=lD,ncol=n)
  wc_hat_s_m <- rep(0,n)
  wctmp <- hat_s_m
  for (i in 1:length(D)){
    wc_hat_s_m[1:(2^i)] <- wc[1:(2^i)]-res$sigma_v*2^(-i/2)
    hat_s_m[i,] <- IWT_PO(wc_hat_s_m,j0,qmf)
    wctmp[i,] <- wc_hat_s_m

    wc_odd_lin[1:(2^i)]  <- wc_odd[1:(2^i)]
    wc_even_lin[1:(2^i)] <- wc_even[1:(2^i)]
    hat_s_odd[i,]  <- IWT_PO(wc_odd_lin,j0,qmf)
    hat_s_even[i,] <- IWT_PO(wc_even_lin,j0,qmf)
  }
  bar_s_odd <- 0.5*(hat_s_odd[,1:(n1-1)] + hat_s_odd[,2:n1])
  bar_s_odd <- cbind(bar_s_odd,hat_s_odd[,1])
  bar_s_even <- 0.5*(hat_s_even[,1:(n1-1)] + hat_s_even[,2:n1])
  bar_s_even <- cbind(bar_s_even,hat_s_even[,1])

  mat_Y_odd <- repmat(Y_odd,length(D),1)
  mat_Y_even <- repmat(Y_even,length(D),1)

  tmp <- (mat_Y_odd-bar_s_even)^2 +
         (mat_Y_even-bar_s_odd)^2
  Crit_CV <- rowSums( tmp )
  return(list("CritCV"=Crit_CV,"hat_s_m"=hat_s_m,"wctmp"=wctmp))
}

