## Collection of GenSequences in a classification tree
## Creator function creates an empty collection

GenCollection <- function(classification=GenClass16S(), 
	type=c("NSV", "sequence"), annotation=NA) 
{
    type <- match.arg(type)
    #x is the collection object, data contains sequences, 
    x <- list(data=list(),classification=list(classification), 
	    type= type, annotation=annotation)      
    class(x) <- "GenCollection"
    return(x)
    
}

list.GenCollection <- function(object, location)
{
    ### FIXME: like count but without the counting

}


count.GenCollection <- function(object, location)
{
    
    levels <- names(object$classification)
    m <- pmatch(level,levels)
    if (is.na(m)) stop("error in level")
    
    ## FIXME: Tree
    l <-array()    
    for(i in 1:length(object))
    {
      l[i]<-object[[i]]$classification[m]
      
    }
        
    t<-table(l)
    #print(t)
    #l<- unique(l)
    return(t)
}


select.GenCollection <- function(x, location) {
    stop("Not implemented!")
}


getGenSequences.GenCollection(x) {
}


print.GenCollection <- function(object) {
    
    cat("Object of class GenCollection\n")
    
    ## report some basic information
}


## gets a GenCollection of type sequence and returns the same structure 
## with NSV (counts)
toNSV.GenCollection <- function(x) {
    stop("Not implemented")
}
