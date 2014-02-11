GenModel <- function(x, rank=NULL, name = NULL, 
	measure="Manhattan", threshold=30, 
	saveClusterInfo=FALSE) {

    d <- .make_stream(x)
    if(measure=="Kullback") d <- d + 1

    emm <- EMM(measure=measure,threshold=threshold)
    build(emm, d)
    
    ### cluster info
    if (saveClusterInfo) clusterInfo <- .getClusterInfo(list(), 
	    last_clustering(emm), attr(x,"id"), 0)
    else clusterInfo <- NA

    structure(list(name=name, rank=rank, nSequences=length(x),
	    model=emm, measure=measure, threshold=threshold, 
	    window=attr(x,"window"), 
	    word=attr(x,"word"), 
	    overlap=attr(x,"overlap"),
	    last_window=attr(x,"last_window"),
	    clusterInfo=clusterInfo), class="GenModel")
}

# creates an model from sequences in the db 
GenModelDB <- function(db, rank=NULL, name=NULL, table="NSV", 
                       measure="Manhattan", threshold=30, 
                       selection=NULL, limit=NULL, random=FALSE, 
		       saveClusterInfo=FALSE, prune=1) {
  
  rank <- BioTools::.pmatchRank(db, rank)
  
  if(random && !is.null(selection)) 
    stop("Cannot use random and selection!")
  
  #check if table exists in db
  if (!(table %in% listGenDB(db)))
    stop("Could not find table in database")
  
  if(!is.null(selection))  nSequences <- length(selection)
  else{
    nSequences <- nSequences(db, rank, name, table=table)
    if (!is.null(limit)) nSequences <- min(nSequences,limit)  
  }
    
  #check for random (selection is NULL)
  if(random) {   
    ids <- getIDs(db, whereRank=rank, whereName=name)
    selection <- sample(ids, nSequences)
  }
  
  
  #emm
  emm <- EMM(measure=measure,threshold=threshold)
  

  #clusterinfo stores the last_clustering details
  if(saveClusterInfo) clusterInfo <- list(nSequences)
  else clusterInfo <- NA
  
  # loop
  i<-0
  total<-0
  while(i<nSequences){
    
    #get 100 sequences at a time
    toGet <- min(100, nSequences-i) 
    if(is.null(selection))
      d <- getSequences(db, rank, name, table=table, limit=c(i,toGet))
    else	
      d <- getSequences(db, rank="id", name=selection, 
                        table=table, limit=c(i,toGet))
    
    
    ids <- names(d)
    d <- .make_stream(d)
    
    # Kullback can not handle 0 counts!
    if(measure=="Kullback") d <- d + 1
    
    #get actual number of sequences
    build(emm, d)

    if(saveClusterInfo) 
      clusterInfo <- .getClusterInfo(clusterInfo,
                                     last_clustering(emm), ids, i)

    
    if(prune>0) prune(emm, count=prune, copy=FALSE, compact=FALSE)
    
    reset(emm) ### make sure there is a NA here
    
    #update value of i 
    i <- i + length(ids)
    cat("GenModel: Processed",i,"sequences\n")
  }
  
  
  if(!is.null(name) && !is.null(rank)) {
    hierarchy <- getHierarchy(db, rank, name)
    if (length(name) > 1)
      hierarchy <- as.data.frame(hierarchy)
    rankName <- hierarchy[rank]
  }else{
    rankName <- NA
    hierarchy <- NA
  }
  
  #get metadata about the table
  meta <- as.character(subset(metaGenDB(db),name==table)["annotation"])
  meta <- sapply(unlist(strsplit(meta,";")), strsplit, "=")  
  meta <- structure(sapply(meta, "[", 2), names=sapply(meta, "[", 1)) 
  
  genModel <- list(name=rankName, rank=rank, 
                   nSequences=nSequences, model=emm, 
                   hierarchy=hierarchy, measure=measure, threshold=threshold, 
                   window=as.integer(meta["window"]),
                   overlap=as.integer(meta["overlap"]), 
                   word=as.integer(meta["word"]), 
                   last_window=as.logical(meta["last_window"]))

  if (saveClusterInfo) genModel$clusterInfo <- clusterInfo
  class(genModel) <- "GenModel"	
  
  genModel		
}

#	Returns the model states and the ID and the segment number of the sequences that are part of that state in format id:segment.
#	By default, returns all states as a list, if a modelState is specified. returns only the sequences that are part of that state
getModelDetails <- function(model, state=NULL, db=NULL)
{	
    if(is.null(model$clusterInfo) || is.na(model$clusterInfo))
	stop("Model needs to be built with saveClusterInfo=TRUE!")
    
    rank = model$rank 
    if(!is.null(state)) {
	occ <- lapply(model$clusterInfo, FUN=function(y) which(y==state))
	occ[sapply(occ, length) <1] <- NULL
	if (!is.null(db)) {     
	    return(
		    data.frame(sequence=	
			    rep(names(occ), times=sapply(occ, length)),
			    segment=unlist(occ), rank=sapply(rep(names(occ), times=sapply(occ, length)), 
                                           FUN=function(x) getRank(db,rank=rank, 
                                          whereRank="id",whereName=x, all=TRUE)),
			    row.names=NULL, stringsAsFactors=FALSE)
		    )
	}
	else {
	    return(
		    data.frame(sequence=	
			    rep(names(occ), times=sapply(occ, length)),
			    segment=unlist(occ), 
			    row.names=NULL, stringsAsFactors=FALSE)
		    )

	}

    }

    l <- lapply(clusters(model$model), FUN=function(x) 
	    getModelDetails(model=model, state=x, db=db))

    names(l) <- clusters(model$model)
    l
}


#	Returns the sequences that are part of a given model state as list of DNA or NSV sequence objects
#
getModelSequences <- function(db, model, state, table="sequences")
{	
    if(is.null(model$clusterInfo) || is.na(model$clusterInfo))
	stop("Model needs to be built with saveClusterInfo=TRUE!")
    
    #get the ids that are part of the model state as a list
    ids<-getModelDetails(model, state, db)

    sequence <- ids[,"sequence"]
    segment <- as.integer(as.character(ids[,"segment"]))
    window <- as.numeric(model$window)
    start <- (segment - 1) * window + 1

    stateSequences <- lapply(1:length(sequence), 
	    FUN=function(i) getSequences(db, rank="id",
		    name=sequence[i], 
		    table, start = start[i], length=window))
    
    do.call(c, stateSequences)
}

#takes a model and returns the clusterInfo i.e. how many states and which sequence goes to which cluster
.getClusterInfo <- function(clusterInfo, last_clustering, ids, offset)
{
    #states are separated by NAs, so get position of NAs
    startStates <- which(is.na(last_clustering))
    for(j in 1:length(startStates))
    {
	start<-startStates[j]+1
	#check if its the last one
	if(is.na(startStates[j+1])) 
	    end<-length(last_clustering)
	else
	    end<-startStates[j+1]-1
	
	clusterInfo[[offset+j]]<-as.integer(last_clustering[start:end])
    }
   
    names(clusterInfo)[(offset+1):(offset+length(ids))] <- ids
    
    return(clusterInfo)
}

### print basic info about a model
print.GenModel <- function(x, ...) {
    cat("Object of class GenModel with", x$nSequences, "sequences\n")
    if(!is.null(x$rank) && !is.null(x$name)) {
	cat(x$rank,":", paste(x$name, collapse=",") , "\n", sep=" ")
    }
    cat("\nModel:\n")
    print(x$model)
}

### plot a model
plot.GenModel <- function(x, ...) {
    #plot(x$model, main=paste(x$rank,": ", x$name, sep=""), ...)
    ### main does not work for interactive!
    plot(x$model, ...)
}


### prune is a S4 generic 
setOldClass("GenModel")
setMethod("prune", signature(x = "GenModel"),
	function(x, ...) {
	    x$model <- prune(x$model, ...)
	    x
	})

recluster <- function(x, method=recluster_kmeans, ...) {
    new_model <- method(x$model, ...)
    x$model <- new_model
    ### FIXME: Is there something to fix in x
    x
}

