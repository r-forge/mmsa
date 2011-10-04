
GenModel <- function(x, threshold=30, measure="manhattan", ...) {
    data <- NSV_make_stream(x$sequences)
    model <- EMM(threshold=threshold, measure=measure, data=data)
    sequences <- matrix(rep(x$class, length(x)), 
	    ncol=length(x$class), byrow=TRUE, 
	    dimnames=list(NULL, names(x$class)))

    l <- list(model=model, sequences=sequences, ...)
    class(l) <- "GenModel"
    l
}


plot.GenModel <- function(x, y=NULL, method="igraph", ...) {
    plot(x$model, method=method)
}
