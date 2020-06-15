# fig2_bot - produces Figure 2. bottom panel of the paper:
#
# Nonparametric estimation in a regression model with additive
# and multiplicative noise.
# Journal of Computational and Applied Mathematics, Vol. 380
#
#Copyright (c) 2020 C. Chesneau, S. El Kolei, J. Kou and F. Navarro

rm(list=ls())
library(rwavelet)
source("makedata.R")
source("rescale.R")
source("twoFCVlinear.R")
source("twoFCVnonlinear.R")
source("PlotWaveCoeffScale.R")

Pn <- function(x0,x){apply((x-x0)^2,1,mean)}
mise <- function(x,y){
  d=mean((as.vector(x)-as.vector(y))^2)
  return(d)
}

signal_names <- c("Blip")
n <- 2^12
J <- log2(n)
j0 <- 0
D  <- 2^(1:(J-1))
filter <- 'Daubechies'
qmf <- MakeONFilter(filter,8)

for (i in 1:length(signal_names)) {
  set.seed(0)
  res <- makedata(n,signal_names[i])
  Y <- res$Y^2/res$sigma_u
  f <- res$f^2
  # noisy wavelet coefficients
  wc <- FWT_PO(Y,j0,qmf)
  # original wavelet coefficients
  wcf <- FWT_PO(f,j0,qmf)
  # 2FCV linear
  cv <- twoFCVlinear(n,Y,j0,qmf,D,wc,res$sigma_v)
  hat_m_2FCV <- which.min(cv$CritCV)
  hat_s_m_2FCV <- cv$hat_s_m[hat_m_2FCV,]
  # Oracle linear
  err <- Pn(cv$hat_s_m,repmat(f,length(D),1))
  m_star <- which.min(err)
  hat_s_m_star <- cv$hat_s_m[m_star,]
  # compute wavelet coefficient using hat_m_2FCV as coarset scale
  wct <- FWT_PO(Y,hat_m_2FCV,qmf)
  wco <- FWT_PO(Y,m_star,qmf)
  # grid value for threshold parameter
  lambda <- 2*sqrt(2*log(n))*0.01# res$sigma_v
  threshval <- seq(lambda/10,100*lambda,length.out=100)
  scale <- (m_star-1):(J-1)
  wch <- wct
  wcc <- array(0,dim=c(length(threshval),length(scale)))
  for (il in 1:length(scale)){
    for (ii in 1:length(threshval)){
      wch[dyad(scale[il])] <- HardThresh(wct[dyad(scale[il])],
                                         threshval[ii])
      wcc[ii,il] <-  mise(wcf,wch)
    }
  }
  indthresh <- apply(wcc,2,which.min)
  wchh <- wct
  for (il in 1:length(scale)){
    wchh[dyad(scale[il])] <- HardThresh(wct[dyad(scale[il])],
                                        threshval[indthresh][il])
  }
  hat_f_oracle_lev <- IWT_PO(wchh,hat_m_2FCV,qmf)
  # Global Hard thresholding
  j0 <- hat_m_2FCV
  wct <- FWT_PO(Y,j0,qmf)
  test <- wct
  test[(2^(j0)+1):n] <- HardThresh(wct[(2^(j0)+1):n],3)
  hat_test <- IWT_PO(test,j0,qmf)
  thresh <- twoFCVnonlinear(n,Y,j0,qmf,wct,threshval)
  
  hat_lambda_2FCV <- which.min(thresh$CritCV)
  hat_s_hard_m_2FCV <- thresh$hatfH[hat_lambda_2FCV,]
}

PlotWaveCoeffScale(wcf,0,1)
PlotWaveCoeffScale(wc,0,1)
PlotWaveCoeffScale(cv$wctmp[hat_m_2FCV,],0,1)
