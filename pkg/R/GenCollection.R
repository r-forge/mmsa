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

findNSV<- function(x)
{
  
  finalNSV <- NULL
  if(length(names(x))==0)
  {
    return(x)
  }
    
  else
  {
    for(i in 1:length(x))
    {
      tempNSV<- findNSV(x[[i]])
      #cat("location = ",location)
      if (length(tempNSV)>0)
      {
        finalNSV<- c(finalNSV,tempNSV)  
      }
      
    }
    #return(c(sapply(x, FUN= function(y) findNSV(y, location-1))))    
  }
  return(finalNSV)
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
toNSV.GenCollection <- function(x) {
    stop("Not implemented")
}
