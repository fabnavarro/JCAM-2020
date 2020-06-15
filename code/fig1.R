# fig1 - produces Figure 1. of the paper:
#
# Nonparametric estimation in a regression model with additive
# and multiplicative noise.
# Journal of Computational and Applied Mathematics, Vol. 380
#
#Copyright (c) 2020 C. Chesneau, S. El Kolei, J. Kou and F. Navarro

rm(list=ls())
source("makedata.R")
source("rescale.R")

signal_names <- c("Parabolas", "Ramp", "Blip")
n <- 4096
for (i in 1:length(signal_names)) {
  set.seed(0)
  res <- makedata(n,signal_names[i])
  plot(res$X, res$f^2, type="l", col = "black", lwd=2,
       axes=FALSE, ann=FALSE, lty=1)
  axis(1, at=seq(0,1,length.out = 6))
  axis(2, las=0)
  box()
  title(xlab="X",cex.lab=1.2, col.lab=rgb(0,0,0))
  title(ylab="r",cex.lab=1.2, col.lab=rgb(0,0,0))
}
