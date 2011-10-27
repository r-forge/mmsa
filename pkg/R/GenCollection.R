#source dependencies

source("Greengenes.R")
## Collection of GenSequences in a classification tree
## Creator function creates an empty collection

GenCollection <- function(classification=GenClass16S_Greengenes(), 
	type=c("NSV", "sequence"), annotation=NA) 
{
    type <- match.arg(type)
    #x is the collection object, data contains sequences, 
    #data refers to sequences or NSV data?
    x <- list(data=list(),classification=list(classification), 
	    type= type, annotation=annotation)      
    class(x) <- "GenCollection"
    return(x)
    
}

list.GenCollection <- function(object, location)
{
    ### FIXME: like count but without the counting
  
    levels <- c(as.vector(names(object$classification[[1]])),"org_name")
    index <- which(levels==location,arr.ind=T)
    #if(length(index)<0)
    #  stop("error in level")
    
    m <- pmatch(location,levels)
    if (is.na(m)) stop("error in level")
    
   totalElements <- length(object$classification)
   if (length(names(object$data[[1]]))==totalElements)
        return(totalElements)
    else
      #recursion
    {
      
        nameList<-findList(object$data,m)
        cat(location," Levels are :\n")
        for(k in 1:length(nameList))
          cat(k," ",nameList[k],"\n")
        
      
    }

}

#helper function
findList <- function(x,location,currLocation=1)
{
  if (currLocation ==location)
  {
      listLast<-NULL
      #for(j in 1:length(x))
      #{
      #  listLast<-c(listLast,names(x[j]))        
      #}
      
      #return (listLast)
      return(names(x))
      }
  else
    #recursion
  {
    currLocation <- currLocation + 1
    nameList<-NULL
    for(i in 1:length(x))    
    {        
      tmpList<-findList(x[[i]],location,currLocation)
      nameList = c(nameList,tmpList)      
    }
  }
  return(nameList)
}


count.GenCollection <- function(object, location)
{
    
    levels <- c(as.vector(names(object$classification[[1]])),"org_name")
    index <- which(levels==location,arr.ind=T)
    
    m <- pmatch(location,levels)
    if (is.na(m)) stop("error in level")
    
   totalElements <- length(object$classification)
   if (length(names(object$data[[1]]))==totalElements)
        return(totalElements)
    else
      #recursion
    {
  
        count=findCount(object$data,m)
        cat("There are  ",count," levels\n")
      
    }
}

#helper function
findCount2 <- function(x,location,currLocation=1)
{
  count=0
  if (currLocation ==location)
  {
      #countLast=0
      #for(j in 1:length(x))
      #{
      #  countLast=cophyuntLast+1
      #}
      #return (countLast)
      return(length(x))
  }
  else
    #recursion
  {
    currLocation <- currLocation + 1
    for(i in 1:length(x))    
    {        
      tmpCount<-findCount(x[[i]],location,currLocation)
      count = count + tmpCount
    }
  }
  return(count)
}

findCount <- function(x,location)
{
  if (location<=0) return(length(x))
  else {
    #recursion
    return(sum(sapply(x, FUN= function(y) findCount(y, location-1))))
  
  }
}


select.GenCollection <- function(x, location) {
    stop("Not implemented!")
}


getGenSequences.GenCollection <- function(x) {
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
