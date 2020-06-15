PlotSpikesScale <- function(base, t, x, L, J) {
  tt <- rbind(t, t, t)
  b <- rep(0, length(x)) + base
  xx <- rbind(b, x + base, b)
  u <- cbind(0, as.vector(tt), 1)
  v <- cbind(base, as.vector(xx), base)
  return(plot(u, v, type = "l", xlim = c(0, 1), ylim = c(-J, -L + 1), axes = FALSE,
              xlab = "", ylab = "", lwd=2))
}
