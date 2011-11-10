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


.findLevel <- function(x, rank) {
    if (is.numeric(rank)) return(rank)

    m <- pmatch(tolower(rank),tolower(names(x$classification)))
    if (is.na(m)) stop("Error in level")
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
	else {
	    #recursion
	    return(unlist((sapply(x, FUN= function(y) findList(y, rank-1)))))
	}
    }


    m <- .findLevel(object, rank)

    #recursion
    nameList<-.findList(object$data,m)
    names(nameList)<- seq(1:length(nameList))
    return(nameList)        
}


# number of nodes at rank
nNodesLevel <- function(object, rank) length(showLevel(object, rank))

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


#helper function

## finds a subtree
## location is a numeric vector with the indices for the ranks
## FIXME: translation from names to location
.getSubTree <- function(x, location) x[[location]]

## returns a location
.findSubTreeLocation <- function(x, rank, name)
{
	.rec <- function(x, rank)
	if (rank<=1) {
	    pmatch(name,names(x))
	}
	    else {
	    #recursion
	    ret <- lapply(x, FUN= function(y) .rec(y, rank-1))
	    found <- sapply(ret, FUN = function(y) !is.na(y[1]))
	    
	    if(any(found)) {
		## FIXME: check for multiple matches!
		return(c(which(found), ret[[which(found)]]))
	    } else return(NA)
	}

	.rec(x, rank)
}




#returns a subtree
findSubTree <- function(x, rank)
{
  levelList <- NULL
  if (rank<=1) {
      return(x)
    } 
    
  else {
    #recursion      
       for(i in 1:length(x))
       {
         ll<-findSubTree(x[[i]],rank-1)
         levelList <- c(levelList,ll)
       }
  }
  return(levelList)
}

.getHeight <- function(x) {
    height <- 0
    while(is(x, "list")) {
	x <- x[[1]]
	height <- height+1
    }   
    height
}

findLeaves<- function(x) .unlist(x, 1, .getHeight(x))    



findLeavesNSV__remove<- function(x)
{
  
  finalLeaves<-list()
  count <-0
   if(names(x)[[1]]=="sequences")
  {
	print(class(x$sequences))   
	return(x$sequences)
  }
    
  else
  {
    for(i in 1:length(x))
    {
      tempLeaf<- findLeaves(x[[i]])
      if (length(tempLeaf)>0)
      {
	finalLeaves<-c(finalLeaves,tempLeaf)
      }
      
    }
        
  }
  return(finalLeaves)
}




#select.GenCollection <- function(x, rank) {
#      d <- x$data[[rank]]  
#      d2 <- list()
#      for (i in 1:length(rank)) {
#        d2 <- list(d2)
#      }
#      d2[[rep(1, length(rank))]] <- d
#      
#      x$data <- d2
#      x
#}
#

getGenSequences.GenCollection <- function(x) {
    .unlist(x, 1, length(x$classification))
  
}

.unlist <- function(x, level, maxLevel) {
    if(level >= (maxLevel)) return(x)
    
    return(unlist(lapply(x, .unlist, level+1, maxLevel), recursive=FALSE))
}


print.GenCollection <- function(object) {
    
    cat("Object of class GenCollection\n")
    
    ## report some basic information
}


## gets a GenCollection of type sequence and returns the same structure 
## with NSV (counts)
toNSV.GenCollection <- function(object, window=100, overlap=0, last_window=FALSE, word=3) {
    x<- object$data

    .rec <- function(x) {
	if(is(x,"GenSequences"))
	{
	    return(toNSV(x, window=window, overlap=overlap, last_window=last_window, word=word))
	}

	else
	{
	    for(i in 1:length(x))
	    {
		x[[i]]<-.rec(x[[i]])

	    }
	    return(x)
	}
    }

    object$data <- x
    object$type <- "NSV"

    object
}


