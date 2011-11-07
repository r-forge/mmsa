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

GenCollectionNew <- function()
{
    x<- list()
    class(x) <- "GenCollection"
    return(x)
}

list.GenCollection <- function(object, location)
{
     
    levels <- as.vector(names(object$classification[[1]]))
    index <- which(levels==location,arr.ind=T)
       
    m <- pmatch(location,levels)
    if (is.na(m)) stop("error in level")
    

      #recursion
    {
      
        nameList<-findList(object$data,m)
        names(nameList)<- seq(1:length(nameList))
        #cat(location," Levels are :\n")
        #for(k in 1:length(nameList))
        #  cat(k," ",nameList[[k]],"\n")
        return(nameList)        
      
    }

}
#helper function
findList <- function(x, location)
{
  if (location<=1) {
      return(names(x))
    } 
    
  else {
    #recursion
    
    return(unlist((sapply(x, FUN= function(y) findList(y, location-1)))))
    
  }
}




count.GenCollection <- function(object, location)
{
    
    levels <- c(as.vector(names(object$classification[[1]])),"org_name")
    index <- which(levels==location,arr.ind=T)
    
    m <- pmatch(location,levels)
    if (is.na(m)) stop("error in level")
    

      #recursion
    {
  
        count=findCount(object$data,m)
        cat("There are  ",count," levels \n")
      
    }
}

#helper function
findCount <- function(x,location)
{
  if (location<=1) return(length(x)) 
  else {
    #recursion
    return(sum(sapply(x, FUN= function(y) findCount(y, location-1))))  
  }
}

genModel.GenCollection<- function(object,location)
{
  
  levels <- as.vector(names(object$classification[[1]]))
  leafLocation <- length(levels)
  index <- which(levels==location,arr.ind=T)
  m <- pmatch(location,levels)
  if (is.na(m)) stop("error in level")
  #start model
  #findSubTree finds sub-tree at the specified level eg: kingdom, phylum, etc
  locationObjects <- findSubTree(object$data,m)
  for(obj in locationObjects)
  {
    emm <- EMM("Manhattan", threshold=0.10)
    for(i in 1:length(obj))
    {
    #find NSV finds the NSVs for that object
    NSVData <-make_stream(findNSV(obj[[i]]))     
    emm<-build(emm,NSVData)
    plotname<-paste(names(obj[[i]]),i,".pdf",sep="")
    #plot to file
    pdf(plotname)
    plot(emm,main=paste(names(obj[[i]]),i))
    dev.off()
    
    reset(emm)
    }
    
    #print(obj)
  }
  #end model
  
}


#helper function
#returns a subtree
findSubTree <- function(x, location)
{
  levelList <- NULL
  if (location<=1) {
      return(x)
    } 
    
  else {
    #recursion      
       for(i in 1:length(x))
       {
         ll<-findSubTree(x[[i]],location-1)
         levelList <- c(levelList,ll)
       }
  }
  return(levelList)
}

findLeaves<- function(x)
{
  
  finalLeaves <- NULL

   if(names(x)[[1]]=="sequences")
  {
	#cat("Type of Last is ",class(x),"\n")
	#print(x$classification) 
	return(x$sequences)
  }
    
  else
  {
    for(i in 1:length(x))
    {
      tempLeaf<- findLeaves(x[[i]])
      #cat("location = ",location)
      if (length(tempLeaf)>0)
      {
        finalLeaves<- c(finalLeaves,tempLeaf)  	
      }
      
    }
    #return(c(sapply(x, FUN= function(y) findNSV(y, location-1))))    
  }
  return(finalLeaves)
}

findLeavesNSV<- function(x)
{
  
  finalLeaves<-list()
  count <-0
   if(names(x)[[1]]=="sequences")
  {
	#cat("Type of Last is ",class(x),"\n")
	#print(x$classification) 
	print(class(x$sequences))   
	return(x$sequences)
  }
    
  else
  {
    for(i in 1:length(x))
    {
      tempLeaf<- findLeaves(x[[i]])
      #cat("location = ",location)
      if (length(tempLeaf)>0)
      {
	finalLeaves<-c(finalLeaves,tempLeaf)
      }
      
    }
        
  }
  return(finalLeaves)
}




select.GenCollection <- function(x, location) {
      d <- x$data[[location]]  
      d2 <- list()
      for (i in 1:length(location)) {
        d2 <- list(d2)
      }
      d2[[rep(1, length(location))]] <- d
      
      x$data <- d2
      x
}


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
toNSV.GenCollection <- function(x,window=100, overlap=0, last_window=FALSE, word=3) {
    	  #x<- object$data
	  if(names(x)[[1]]=="sequences")
	  {
		cnt <- count_sequences(x$sequences,window=window, overlap=overlap, word=word,last_window=last_window)
      		stream <- make_stream(cnt)
 		x$sequences<-cnt		
		#names(x[[1]])<-"sequences"
		return(x)
	  }
	    
	  else
	  {
	    for(i in 1:length(x))
	    {
	      #tempLeaf<- findLeaves(x[[i]])
	      #cat("location = ",location)
	      #if (length(tempLeaf)>0)
	      #{
		#finalLeaves<- c(finalLeaves,tempLeaf)  
	      #}
	      x[[i]]<-toNSV.GenCollection(x[[i]])
	      
	    }
	    #return(c(sapply(x, FUN= function(y) findNSV(y, location-1))))    
	  }
	  #object$data <- x
	  return(x)
}
