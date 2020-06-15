# fig3 - produces Figure 3. of the paper:
#
# Nonparametric estimation in a regression model with additive
# and multiplicative noise.
# Journal of Computational and Applied Mathematics, Vol. 380
#
#Copyright (c) 2020 C. Chesneau, S. El Kolei, J. Kou and F. Navarro

rm(list=ls())
library(rwavelet)
library(ggplot2)
source("makedata.R")
source("rescale.R")
source("twoFCVlinear.R")
source("twoFCVnonlinear.R")
source("SUREsoft.R")
source("PlotWaveCoeffScale.R")

signal_names <- c("Parabolas", "Ramp","Blip")

Pn <- function(x0,x){apply((x-x0)^2,1,mean)}
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

MC <- 100
MSE_lin_hat <- MSE_non_hat <- matrix(0,nrow =length(signal_names),
                                     ncol = MC)
MSE_lin_ora <- MSE_non_ora <- matrix(0,nrow =length(signal_names),
                                     ncol = MC)

thresh_2FCV <- thresh_oracle <- matrix(0,nrow =length(signal_names),
                                       ncol = MC)

for (i in 1:length(signal_names)) {
  set.seed(0)
  for (it in 1:MC){
    res <- makedata(n,signal_names[i])
    Y <- res$Y^2/res$sigma_u
    f <- res$f^2
    wc <- FWT_PO(Y,j0,qmf)
    wcf <- FWT_PO(f,j0,qmf)
    # linear CV
    cv <- twoFCVlinear(n,Y,j0,qmf,D,wc,res$sigma_v)
    hat_m_2FCV <- which.min(cv$CritCV)
    hat_s_m_2FCV <- cv$hat_s_m[hat_m_2FCV,]
    # linear oracle
    err <- Pn(cv$hat_s_m,repmat(f,length(D),1))
    m_star <- which.min(err)
    hat_s_m_star <- cv$hat_s_m[m_star,]
    
    lambda <- 2*sqrt(2*log(n))*0.01# res$sigma_v
    threshval <- seq(lambda/10,100*lambda,length.out=50)
    # global Hard thresholding
    j0 <- hat_m_2FCV
    wct <- FWT_PO(Y,j0,qmf)
    thresh <- twoFCVnonlinear(n,Y,j0,qmf,wct,threshval)
    hat_lambda_2FCV <- which.min(thresh$CritCV)
    hat_s_hard_m_2FCV <- thresh$hatfH[hat_lambda_2FCV,]
    
    errhard <- Pn(thresh$hatfH,repmat(f,length(threshval),1))
    hat_oracle <- which.min(errhard)
    hat_s_hard_m_oracle <- thresh$hatfH[hat_oracle,]
    
    thresh_2FCV[i,it] <- hat_lambda_2FCV
    thresh_oracle[i,it] <- hat_oracle
    MSE_lin_hat[i,it] <- mise(hat_s_m_2FCV,f)
    MSE_lin_ora[i,it] <- mise(hat_s_m_star,f)
    MSE_non_hat[i,it] <- mise(hat_s_hard_m_2FCV,f)
    MSE_non_ora[i,it] <- mise(hat_s_hard_m_oracle,f)
  }
}
df <- data.frame(t(MSE_lin_hat),
                 t(MSE_lin_ora),
                 t(MSE_non_hat),
                 t(MSE_non_ora))

df1 <- df[c(1,4,7,10)]
names(df1) <- c("2FCVlin","MSElin","2FCVnon","MSEnon")
df2 <- df[c(2,5,8,11)]
names(df2) <- c("2FCVlin","MSElin","2FCVnon","MSEnon")
df3 <- df[c(3,6,9,12)]
names(df3) <- c("2FCVlin","MSElin","2FCVnon","MSEnon")

tmp1 <- data.frame("methods"=rep(names(df1),each=MC),
                   "MSE"=c(df1$`2FCVlin`,
                           df1$MSElin,
                           df1$`2FCVnon`,
                           df1$MSEnon))
ylim1 = boxplot.stats(tmp1$MSE)$stats[c(1, 5)]

tmp2 <- data.frame("methods"=rep(names(df2),each=MC),
                   "MSE"=c(df2$`2FCVlin`,
                           df2$MSElin,
                           df2$`2FCVnon`,
                           df2$MSEnon))
ylim2 = boxplot.stats(tmp2$MSE)$stats[c(1, 5)]

tmp3 <- data.frame("methods"=rep(names(df3),each=MC),
                   "MSE"=c(df3$`2FCVlin`,
                           df3$MSElin,
                           df3$`2FCVnon`,
                           df3$MSEnon))
ylim3 = boxplot.stats(tmp3$MSE)$stats[c(1, 5)]

ggplot(data = tmp1, aes(x = methods, y = MSE, fill = methods))+
  stat_boxplot(geom ='errorbar', width = 0.5)+
  geom_boxplot(stat = "boxplot",outlier.colour="red")+
  theme(legend.position="none",
        axis.text=element_text(size=13),
        axis.title.y=element_text(size=14),
        axis.title.x=element_blank())

ggplot(data = tmp2, aes(x = methods, y = MSE, fill = methods))+
  stat_boxplot(geom ='errorbar', width = 0.5)+
  geom_boxplot(stat = "boxplot",outlier.colour="red")+
  theme(legend.position="none",
        axis.text=element_text(size=13),
        axis.title.y=element_text(size=14),
        axis.title.x=element_blank())

ggplot(data = tmp3, aes(x = methods, y = MSE, fill = methods))+
  stat_boxplot(geom ='errorbar', width = 0.5)+
  geom_boxplot(stat = "boxplot",outlier.colour="red")+
  theme(legend.position="none",
        axis.text=element_text(size=13),
        axis.title.y=element_text(size=14),
        axis.title.x=element_blank())
