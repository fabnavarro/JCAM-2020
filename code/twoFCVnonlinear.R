twoFCVnonlinear <- function(n,Y,j0,qmf,wc,threshval){
  wchard <- wc
  lenthresh <- length(threshval)
  n1 <- n/2
  wc_odd_lin  <- rep(0,n1)
  wc_even_lin <- rep(0,n1)
  Y_odd <- Y[seq(1,n,2)]
  Y_even <- Y[seq(2,n,2)]

  wc_odd  <- FWT_PO(Y_odd,j0,qmf)
  wc_even <- FWT_PO(Y_even,j0,qmf)
  hat_s_odd  <- matrix(0,nrow=lenthresh,ncol=n1)
  hat_s_even <- matrix(0,nrow=lenthresh,ncol=n1)

  hatfH <- array(0,dim=c(lenthresh,n))
  wctmp <- hatfH
  for (ii in 1:lenthresh){
    wchard[(2^(j0)+1):n] <- HardThresh(wc[(2^(j0)+1):n],threshval[ii])
    hatfH[ii,] <- IWT_PO(wchard,j0,qmf)
    wctmp[ii,] <- wchard

    wc_odd_lin[(2^(j0)+1):n1]  <- HardThresh(wc_odd[(2^(j0)+1):n1],threshval[ii])
    wc_even_lin[(2^(j0)+1):n1] <- HardThresh(wc_even[(2^(j0)+1):n1],threshval[ii])
    hat_s_odd[ii,]  <- IWT_PO(wc_odd_lin,j0,qmf)
    hat_s_even[ii,] <- IWT_PO(wc_even_lin,j0,qmf)
  }
  # 2-fold CV
  bar_s_odd <- 0.5*(hat_s_odd[,1:(n1-1)] + hat_s_odd[,2:n1])
  bar_s_odd <- cbind(bar_s_odd,hat_s_odd[,1])
  bar_s_even <- 0.5*(hat_s_even[,1:(n1-1)] + hat_s_even[,2:n1])
  bar_s_even <- cbind(bar_s_even,hat_s_even[,1])

  mat_Y_odd <- repmat(Y_odd,lenthresh,1)
  mat_Y_even <- repmat(Y_even,lenthresh,1)

  tmp <- (mat_Y_odd-bar_s_even)^2 +
         (mat_Y_even-bar_s_odd)^2
  Crit_CV <- rowSums(tmp)
  return(list("CritCV"=Crit_CV,"hatfH"=hatfH,"wctmp"=wctmp))
}


