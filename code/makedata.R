makedata <- function(n, name) {
  X <- sort(runif(n))
  if (name == "Blip") {
    sig <- (0.32 + (0.6 * X) + 0.3 * exp(-100 * ((X - 0.3)^2))) * (X <= 0.8) +
      (-0.28 + (0.6 * X) + 0.3 * exp(-100 * ((X - 1.3)^2))) * (X > 0.8)
  }
  if (name == "Ramp") {
    sig <- X - (X >= 0.37)
  }
  if (name == "Parabolas") {
    pos <- c(0.1, 0.2, 0.3, 0.35, 0.37, 0.41, 0.43, 0.5, 0.7, 0.9)
    hgt <- c(-30, 60, -30, 500, -1000, 1000, -500, 7.5, -15, 7.5)
    sig <- 2 * rep(0, length(X))
    for (j in 1:length(pos)) {
      sig <- sig + hgt[j] * ((X - pos[j])^2) * (X > pos[j])
    }
    sig <- sig + 0.8
  }
  a <- -1
  b <- 1
  sigma_v <- sqrt(0.01)
  sigma_u <- (b-a)^2/12+((b+a)/2)^2
  U <- runif(n,a,b)
  V <- rnorm(n)
  Esigma_v2 <- sigma_v^2
  Y <- U*sig + sigma_v*V
  return(list(X=X,Y=Y,f=sig,sigma_u=sigma_u, sigma_v=Esigma_v2))
}
