# fig2_top - produces Figure 2. top panel of the paper:
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

signal_names <- c("Blip")

Pn <- function(x0,x){apply((x-x0)^2,1,mean)}
repmat <- function(a, n, m){
  a <- matrix(a,
              nrow = 1,
              ncol = length(a))
  matrix(1, n, m) %x% a
}
mise <- function(x,y){
  d=mean((as.vector(x)-as.vector(y))^2)
  return(d)
}
n <- 4096
J <- log2(n)
j0 <- 0
D  <- 2^(1:(J-1))
filter <- 'Daubechies'
qmf <- MakeONFilter(filter,8)
for (i in 1:length(signal_names)) {
  set.seed(0)
  hat_s_m <- matrix(rep(0,length(D)*n),nrow=length(D),ncol=n)
  wc_hat_s_m <- rep(0,n)
  res <- makedata(n,signal_names[i])
  Y <- res$Y^2/res$sigma_u
  f <- res$f^2
  wc <- FWT_PO(Y,j0,qmf)
  
  cv <- twoFCVlinear(n,Y,j0,qmf,D,wc,res$sigma_u)
  hat_m_2FCV <- which.min(cv$CritCV)
  hat_s_m_2FCV <- cv$hat_s_m[hat_m_2FCV,]
  
  err <- Pn(cv$hat_s_m,repmat(f,length(D),1))
  m_star <- which.min(err)
  hat_s_m_star <- cv$hat_s_m[m_star,]
  
  #lambda <- 2*sqrt(2*log(n))*res$sigma_v
  #threshval <- seq(lambda/5,50*lambda,length.out=300)
  threshval <- seq(1.5,2.5,length.out=300)
  wct <- FWT_PO(Y,hat_m_2FCV,qmf)
  thresh <- twoFCVnonlinear(n,Y,hat_m_2FCV,qmf,wct,threshval)
  hat_lambda_2FCV <- which.min(thresh$CritCV)
  hat_s_hard_m_2FCV <- thresh$hatfH[hat_lambda_2FCV,]
  print(mise(hat_s_hard_m_2FCV,f))
  print(mise(hat_s_m_2FCV,f))
  psnroutH <- apply(thresh$hatfH,1,mise,y=f)
  hatH <- which.min(psnroutH)
  besthatfH <- thresh$hatfH[hatH,]
  
  limit_sup <- max(f,hat_s_m_2FCV,res$Y^2)
  limit_inf <- min(f,hat_s_m_2FCV,res$Y^2)
  plot(res$X, res$Y^2, col = "grey",
       lwd = 1, axes = FALSE, ann = FALSE, lty = 1,
       ylim = c(limit_inf,limit_sup))
  lines(res$X,res$f^2, type="l", col = "black",lwd=2, lty=1)
  axis(1, at=seq(0,1,length.out = 6))
  axis(2, las=0)
  box()
  matlines(res$X,hat_s_m_2FCV,type="l",lwd=2,lty=1,
           col=c("red"))
  title(xlab="X",cex.lab=1.2, col.lab=rgb(0,0,0))
  title(ylab="r",cex.lab=1.2, col.lab=rgb(0,0,0))

  rCrit_CV = rescale(cv$CritCV,min(err),max(err));
  par(mgp = c(1.5, 0.5, 0), tcl = -0.2,
      mar = .1 + c(2.5,2.5,0,0), oma = c(0,0,0,0))
  plot(1:11,err, type="o", col = "black",lwd=2,
       axes=FALSE, ann=FALSE, lty=1)
  axis(1,at =1:11)
  axis(2, las=0)
  box()
  matlines(1:11,rCrit_CV,type="o",lwd=2,lty=1,pch=22,
           col=c("red"))
  title(xlab="j_*",cex.lab=1.2, col.lab=rgb(0,0,0))
  title(ylab="MSE",
        cex.lab=1.2, col.lab=rgb(0,0,0))

  rCrit_CV = rescale(thresh$CritCV,min(psnroutH),max(psnroutH));
  plot(threshval,psnroutH, type="s",
       col = "black",lwd=2,
       axes=FALSE, ann=FALSE, lty=1)
  axis(1)
  axis(2, las=0)
  box()
  matlines(threshval,rCrit_CV,type="s",lwd=2,lty=1,
           col=c("red"))
  title(xlab="lambda", cex.lab=1.2, col.lab=rgb(0,0,0))
  title(ylab="MSE", cex.lab=1.2, col.lab=rgb(0,0,0))
}
