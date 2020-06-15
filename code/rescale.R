rescale <- function(x,a,b){
	m <- min(as.vector(x))
	M <- max(as.vector(x))
	y <- (b-a)*(x-m)/(M-m)+a
	return(y)
}