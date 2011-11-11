## Collection of GenSequences in a classification tree
## Creator function creates an empty collection

GenCollection <- function(classification=GenClass16S_Greengenes(), 
	type=c("sequence", "NSV"), annotation=NA) 
{
    type <- match.arg(type)
    #x is the collection object, data contains sequences, 
    #data refers to sequences or NSV data?
    x <- list(data=list(),classification=classification, 
	    type= type, annotation=annotation)      
    class(x) <- "GenCollection"
    return(x)
    
}


print.GenCollection <- function(object) {
    cat("Object of class GenCollection")
    cat(" of type", object$type) 
    cat(" with", nSequences(object), "sequences.\n")
}



## find a rank as a number
findRank <- function(x, rank) {
    if (is.numeric(rank)) return(rank)

    m <- pmatch(tolower(rank),tolower(names(x$classification)))
    if (is.na(m)) stop("Error in level")
    
    
    names(m) <- names(x$classification)[m]
    m
}

# the rank is a character string containing a classification level
# list the elements at one level in the tree
showRank <- function(object, rank)
{
    #helper function
    .findList <- function(x, rank)
    {
	if (rank<=1) return(names(x))
	return(unlist((sapply(x, FUN = function(y) 
					.findList(y, rank-1)))))
    }


    m <- findRank(object, rank)

    #recursion
    nameList <- .findList(object$data,m)
    if(!is.null(nameList)) names(nameList) <- seq(1:length(nameList))
    return(nameList)        
}


# number of nodes at rank
nNodesRank <- function(object, rank) length(showRank(object, rank))
nSequences <- function(object) nNodesRank(object, 
	length(object$classification))



## location is a numeric vector with the indices for the ranks
## returns a location
findLocation <- function(x, rank, name)
{
    .rec <- function(x, rank){
	if (rank<=1) {
	    m <- pmatch(name,names(x))
	    names(m) <- names(x)[m]
	    return(m)
	
	} else {
	    #recursion
	    ret <- lapply(x, FUN= function(y) .rec(y, rank-1))
	    found <- sapply(ret, FUN = function(y) !is.na(y[1]))

	    if(any(found)) {
		## FIXME: check for multiple matches!
		m <- which(found)
		names(m) <- names(x)[m]
		return(c(m, ret[[which(found)]]))
	    } else return(NA)
	}
    }
    
    .rec(x$data, rank)
}


## FIXME: translation from names to location
.getSubTree <- function(x, location) x[[location]]

## compute the height of a tree
.getHeight <- function(x) {
    height <- 0
    while(is(x, "list")) {
	x <- x[[1]]
	height <- height+1
    }   
    height
}

## unlists a tree 
.unlist <- function(x, level, maxLevel) {
    if(level >= (maxLevel)) return(x)

    return(unlist(lapply(x, .unlist, level+1, maxLevel), 
		    recursive=FALSE))
}


## return all sequences as a list
getSequences<- function(x, location=NULL) {
    if (is.null(location)) data <- x$data
    else data <- .getSubTree(x$data, location)
    
    .unlist(data, 1, .getHeight(data))    
}

.getSequencesData<- function(x, location=NULL) {
    if (is.null(location)) data <- x
    else data <- .getSubTree(x, location)
    
    .unlist(data, 1, .getHeight(data))    
}

## gets a GenCollection of type sequence and returns the same structure 
## with NSV (counts)
toNSV.GenCollection <- function(object, window=100, 
	overlap=0, word=3, last_window=FALSE) {
    
    new <- object
    new$data <- list()
    new$type <- "NSV"

    # get sequences
    sequences <- getSequences(object)
    
    # count sequences
    sequences <- lapply(sequences, toNSV, window, overlap, 
	    word, last_window)

    # build new tree

    for(s in sequences) {
	
	## recreate branch
	for(i in 1:length(s$classification)) {
	    if(is.null(new$data[[s$classification[1:i]]])) {
		new$data[[s$classification[1:i]]] <- list()
	    }
	}
    
	## insert sequence
	new$data[[s$classification]] <- s
    }

    new
}


#location is a numeric vector indicating subtree
genModel.GenCollection<- function(object,location, measure="Kullback", threshold=0.10)
{
	if (object$type != "NSV") 
		stop("Not in NSV format")
	x<-object$data
	subTree<- x[[location]]
	GenSequences <- .getSequencesData(subTree)
	emm<- EMM(measure=measure,threshold=threshold)
	for(GenSequence in GenSequences)
	{
		sequence<-as.data.frame(GenSequence$sequences)		
		emm<-build(emm,sequence+1)
	}
	#remove last element from location
	location1<-location[1:length(location)-1]
	rank<-names(object$classification)[[3]]
	value<-names(x[[location1]])[location[length(location)]]
	#cat("rank: ",rank," value:",value,"\n")
	op<-paste("Rank : ",rank," Value : ",value)
	names(emm)<-op
	emm    
}

#location is a numeric vector indicating subtree
genPlotModel.GenCollection<- function(object,location, measure="Kullback", threshold=0.10)
{
	if (object$type != "NSV") 
		stop("Not in NSV format")
	x<-object$data
	subTree<- x[[location]]
	GenSequences <- .getSequencesData(subTree)
	emm<- EMM(measure=measure,threshold=threshold)
	for(GenSequence in GenSequences)
	{
		sequence<-as.data.frame(GenSequence$sequences)		
		emm<-build(emm,sequence+1)
	}
	#remove last element from location
	location1<-location[1:length(location)-1]
	rank<-names(object$classification)[[3]]
	value<-names(x[[location1]])[location[length(location)]]
	plotName<-paste("plots/",value,".pdf",sep="")
	title<-paste("Rank : ",rank," Value : ",value)
	pdf(plotName)
	plot(emm,main=title)
	dev.off()	
	fileName<-sub("pdf","RData",plotName)
	save(emm,file=fileName)
	emm
}

plot.GenModel<-function(emm,title)
{
	plot(emm,main=title)
}

