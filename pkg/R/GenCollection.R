
## Collection of GenSequences in a classification tree
## Creator function creates an empty collection
GenCollection <- function(classification, 
	type=c("NSV", "sequence"),
	annotation=NA, ...) {

    type <- match.arg(type)
    x <- list(data=list(), type= type, annotation=annotation, ...)
    
    class(x) <- "GenCollection"
    x
}

print.GenCollection <- function(object) {
    cat("Object of class GenCollection\n")
    ## report some basic information
}
    
## adds a GenSequences object to the collection tree at the correct
## position
add.GenCollection <- function(x, newdata) {
    stop("Not implemented!")
}

## gets a GenCollection of type sequence and returns the same structure 
## with NSV (counts)
toNSV.GenCollection <- function(x) {
    stop("Not implemented")
}
