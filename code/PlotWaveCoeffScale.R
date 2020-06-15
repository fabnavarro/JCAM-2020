source("PlotSpikesScale.R")
PlotWaveCoeffScale <- function(wc, L, scal) {
  wavecoef <- ShapeAsRow(wc)
  n <- dyadlength(wavecoef)$x
  J <- dyadlength(wavecoef)$y
  if (scal == 0) {
    scal <- 1/max(abs(wavecoef[(2^L + 1):n]))
  }
  for (j in seq(J - 1, L, -1)) {
    tj <- (0.5:(2^(j) - 0.5))/2^(j)
    PlotSpikesScale(-j, tj, wavecoef[dyad(j)] * scal, L, J)
    par(new = TRUE)
  }
  axis(2, at = seq(-J + 1, -L, 1), labels = seq(J - 1, L, -1),
       cex.axis=1.5)
  axis(1, at = seq(0, 1, 0.2), labels = seq(0, 1, 0.2),
       cex.axis=1.5)
  title(#ylab="Scales",mgp=c(4,1,0),
        cex.lab=1.2, col.lab=rgb(0,0,0))
  box()
}
