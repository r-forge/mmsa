
## Set of sequences or NSVs for a single object (e.g., strain)
GenSequences <- function(sequences, classification, 
	type=c("NSV", "sequence") ,
	annotation=NA, ...) {

    type <- match.arg(type)

    if(!is(sequences, "list")) sequences <- list(sequences)
    l <- list(sequences=sequences, 
	    classification = classification, 
	    type=type,
	    annotation=annotation,
	    ...)

    class(l) <- "GenSequences"
    l
}

print.GenSequences <- function(object) {
    cat("Object of class GenSequences for", object$type,"\n")
    print(object$classification)
    cat("Number of sequences:", length(object$sequences),"\n")
}

summary.GenSequences <- function(object) {
    cat("Object of class GenSequences for", object$type,"\n")
    print(object$classification)
    cat("Number of sequences:", length(object$sequences),"\n")
    #cat("Fragement sizes:\n")
    #print(summary(sapply(object$sequences, length)))
    cat("Annotation:", object$annotation,"\n\n")
}

length.GenSequences <- function(x) {
    if(is(x$sequences, "matrix")) nrows(x$sequences)
    else length(x$sequences)
    }


## gets a GenSequences of type sequence and returns the same structure 
## with NSV (counts)
## FIXME: this has to be toNSV.GenSequences
toNSV <- function(x, window=100, overlap=0, word=3, 
	offset = 0, last_window=FALSE) {
    
    x$sequences <- 
    NSV_list_counter(x$sequences,
	    window=window, overlap=overlap, word=word, offset=offset, 
	    last_window=last_window)
    x$count_info <- c(window=window, overlap=overlap, word=word, offset=offset,
	    last_window=last_window)
    x$type <- "NSV"
    x
}
