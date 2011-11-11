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



genModel.GenCollection<- function(object,location, method, threshold)
{
    
    # getSubtree (location)
    # findLeaves
    # loop to add to EMM
    # return EMM
    
    levels <- as.vector(names(object$classification[[1]]))
  	m <- pmatch(rank,levels)
  	if (is.na(m)) stop("error in level")
	SubTrees <- findSubTree(object$data,m)
	if (!file.exists("plots")){
    		dir.create(file.path(getwd(), "/plots"))
    	}
	setwd(paste(getwd(),"/plots",sep=""))
	for(i in 1:length(SubTrees))
	{
		leaves <- findLeavesNSV(SubTrees[[i]])
		emm <- EMM("Kullback",threshold=0.10)
		for(leaf in leaves)
			{
				emm<- build(emm,leaf+1)
			}
		plotName<-paste(names(SubTrees)[[i]],".pdf",sep="")
		pdf(plotName)
		plot(emm,main=paste(names(SubTrees)[[i]],"level =",rank))
		dev.off()
		fileName<-sub("pdf","RData",plotName)
		save(emm,file=fileName)
		
		reset(emm)	
	}
  #start
	#leaves<- findLeavesNSV(object)
	#emm<- EMM("Kullback",threshold=0.10)
	#for(leaf in leaves)
	#{#
	#	emm<- build(emm,leaf+1)	
		#reset(emm)	
	#}	
  #end   

  #levels <- as.vector(names(object$classification[[1]]))
  #leafLocation <- length(levels)
  #index <- which(levels==rank,arr.ind=T)
  #m <- pmatch(rank,levels)
  #if (is.na(m)) stop("error in level")
  #start model
  #findSubTree finds sub-tree at the specified level eg: kingdom, phylum, etc
  #rankObjects <- findSubTree(object$data,m)
  #for(obj in rankObjects)
  #{
  #  emm <- EMM("Manhattan", threshold=0.10)
  #  for(i in 1:length(obj))
  #  {
  #  #find NSV finds the NSVs for that object
  #  NSVData <-make_stream(findNSV(obj[[i]]))     
  #  emm<-build(emm,NSVData)
  #  plotname<-paste(names(obj[[i]]),i,".pdf",sep="")
  #  #plot to file
  #  pdf(plotname)
  #  plot(emm,main=paste(names(obj[[i]]),i))
  #  dev.off()
    
  #  reset(emm)
  #  }
    
    #print(obj)
  #}
  #end model
	WD<-getwd()
	WD<-sub("/plots","",WD)	
	setwd(WD)
  
}

