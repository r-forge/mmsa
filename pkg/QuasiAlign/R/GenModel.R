GenModel <- function(x, rank=NULL, name = NULL, 
	measure="Manhattan", threshold=30 
    , showClusterInfo=TRUE) {

    d <- .make_stream(x)
	clusterInfo <- list(length(x))
    if(measure=="Kullback") d <- d + 1

    emm <- EMM(measure=measure,threshold=threshold)
    build(emm, d)
    
    l<-last_clustering(emm)
    
	clusterInfo <- .getClusterInfo(clusterInfo,l,0)
	names(clusterInfo) <- attr(x,"id")
    genModel <- list(name=name, rank=rank, nSequences=length(x), model=emm, measure=measure, threshold=threshold, window=attr(x,"window"), word=attr(x,"word"), overlap=attr(x,"overlap"),
		last_window=attr(x,"last_window"))

    if (showClusterInfo) genModel$clusterInfo <- clusterInfo
    
    class(genModel) <- "GenModel"	
    genModel		
}

# creates an model from sequences in the db 
GenModelDB <- function(db, rank=NULL, name=NULL, table="NSV", 
	measure="Manhattan", threshold=30, 
	selection=NULL, limit=NULL, random=FALSE, showClusterInfo=TRUE) {

    #check if table exists in db
    if (length(which(table==listGenDB(db))) == 0)
	stop("Could not find table in database")

    #check if table is of type NSV	
    meta<-dbReadTable(db$db,"metaData") #meta = table data in memory
    index<-which(meta$name==table)  #find index of table
    #if (meta$type[index]!="NSV")
	#	stop("Not an NSV table")
    #get metadata about the table
    meta<-as.character(subset(metaGenDB(db),name==table)["annotation"])
    x<-unlist(strsplit(meta,";"))	
    window <-as.integer(sub("window=","",x[3]))
    overlap <- as.integer(sub("overlap=","",x[4]))
    word <- as.integer(sub("word=","",x[5]))
    last_window <- as.logical(sub("last_window=","",x[6]))	

    emm <- EMM(measure=measure,threshold=threshold)

    nSequences <- nSequences(db, rank, name, table=table)
    hierarchy <- getHierarchy(db, rank, name)
	if (length(name) > 1)
		hierarchy <- as.data.frame(hierarchy)
	
    if (!is.null(limit)) nSequences <- min(nSequences,limit)
	#check for random and if so get random 'limit' sequences from the DB
	if (random) 
		{	#get the IDs
			ids <- getRank(db, rank="id", whereRank=rank, whereName=name)
			#if(is.null(limit)) limit <- nSequences
			selection <- sample(as.vector(ids),nSequences)
		}
    #clusterinfo stores the last_clustering details
    clusterInfo<-list(nSequences)
    ### Kullback can not handle 0 counts!
    if(measure=="Kullback") d <- d + 1
	
	allRanks <- getRank(db,rank=rank)
	name <- allRanks[pmatch(name,allRanks)]
    cat("GenModel: Creating model for ",rank,": ",paste(name, collapse=",") ,"\n",sep="")
	
    if (is.null(selection)){
		i<-0
		total<-0
		while(i<nSequences){
	    	#get 100 sequences at a time
			d <-getSequences(db,rank,name,table,limit=c(i,100))
	    	n <- length(d)
	    	n <- min(n,100)
	    	ids <- names(d)[1:n]
	    	d <- .make_stream(d)
	    	#get actual number of sequences
	    	build(emm, d)
	    	l<-last_clustering(emm)
	    	clusterInfo<- .getClusterInfo(clusterInfo,l,i)
	    	names(clusterInfo)[(i+1):(i+n)] <- ids
	    	reset(emm)
	    	#update value of i 
	    	i<-min(i+100,nSequences)
	    	cat("GenModel: Processed",i,"sequences\n")
		}
    } else if (!is.null(selection)) {
			d<-getSequences(db, rank="id", name=selection, table, limit=limit)
			ids <- names(d)
			if (length(d)==0) stop("GenModel called with 0 sequences")
			d <- .make_stream(d)
			build(emm, d)
			l<-last_clustering(emm)
			clusterInfo <- .getClusterInfo(clusterInfo,l,0)
			names(clusterInfo) <- ids
			reset(emm)
			nSequences <- length(selection)
			cat("GenModel: Processed ",nSequences ," sequences\n")
    }

    rank <- .pmatchRank(db, rank)
    rankName <- hierarchy[rank]

    genModel <- list(name=rankName, rank=rank, 
	    nSequences=nSequences, model=emm, 
	    hierarchy=hierarchy, measure=measure, threshold=threshold, 
	    window=window, word=word, overlap=overlap, 
	    last_window=last_window)

    if (showClusterInfo) genModel$clusterInfo <- clusterInfo
	rm(db)
    class(genModel) <- "GenModel"	
    genModel		
}

#	Returns the model states and the ID and the segment number of the sequences that are part of that state in format id:segment.
#	By default, returns all states as a list, if a modelState is specified. returns only the sequences that are part of that state
getModelDetails <- function(model, state=NULL, db=NULL)
{	
   	rank = model$rank 
	if(!is.null(state)) {
	occ <- lapply(model$clusterInfo, FUN=function(y) which(y==state))
	occ[sapply(occ, length) <1] <- NULL
	if (!is.null(db)) {     
			return(
				data.frame(sequence=	
					rep(names(occ), times=sapply(occ, length)),
					segment=unlist(occ), rank=sapply(rep(names(occ), times=sapply(occ, length)), FUN=function(x) getRank(db,rank=rank, whereRank="id",whereName=x, all=TRUE,partialMatch=FALSE)),
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
    #get the ids that are part of the model state as a list
    ids<-getModelDetails(model, state, db)

    if (table=="sequences")
		stateSequences <- DNAStringSet()
	else
		stateSequences <- list()
    
    for(i in 1:nrow(ids)){
	id <- ids[i,]
	sequence <- id["sequence"]
	segment <- as.integer(as.character(id["segment"]))
	window <- as.numeric(model$window)
	start <- (segment - 1) * window + 1
	
	stateSequences <- c(stateSequences,
		getSequences(db, rank="id",name=sequence, 
			table, start = start, length=window, partialMatch=FALSE))
    }
    
    #if(table == "sequences") 
	#stateSequences <- do.call(c, stateSequences)
    if (table !="sequences")
		class(stateSequences) <- "NSVSet"

    return(stateSequences)
}

#takes a model and returns the clusterInfo i.e. how many states and which sequence goes to which cluster
.getClusterInfo <- function(clusterInfo, last_clustering, offset)
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



# reads all fasta files in a directory into a db and 
# creates NSV table with all sequences
processSequences <- function(dir, db, metaDataReader=GreengenesMetaDataReader,
	...) {
    for(f in dir(dir, full.names=TRUE))
    {
	cat("Processing file: ",f,"\n")
	addSequences(db, f, metaDataReader=metaDataReader)
    }

    createNSVTable(db, "NSV", ...)
}



# Creates models in modelDir directory for all names in rank.  If selection is
# specified, then it uses only those sequences for creating the model
createModels <- function(modelDir, rank = "Phylum", db, selection=NULL, 
	limit=NULL, ...) 
{

    ### make sure the directory is always lower case
    rank <- tolower(rank)

    ### check if modelDir exists
    ### create rank subdir
    rankDir<-file.path(modelDir,rank)
    if (file.exists(modelDir))
    {
	if(!file.exists(rankDir)) dir.create(rankDir)
    }
    else
    {
	dir.create(modelDir)
	dir.create(rankDir)
    }
    #get All ranks
    rankNames <- getRank(db, rank)
    for(n in rankNames) {
	emm <- GenModelDB(db, table="NSV", rank, name=n,
		selection=selection, limit=limit, ...)
	n<-gsub("/","",n)
	saveRDS(emm, file=paste(rankDir, "/", n, ".rds", sep=''))
    }
}






