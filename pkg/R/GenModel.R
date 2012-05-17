GenModel <- function(x, rank=NULL, name = NULL, 
	measure="Manhattan", threshold=30) {
    
    d <- .make_stream(x)
    if(measure=="Kullback") d <- d + 1

    emm <- EMM(measure=measure,threshold=threshold)
    build(emm, d)
    emm

    genModel <- list(name=name, rank=rank, nSequences=length(x), model=emm)
    
    class(genModel) <- "GenModel"	
    genModel		
}

# creates an model from sequences in the db 
GenModelDB <- function(db, rank=NULL, name=NULL, table="NSV", 
	measure="Manhattan", threshold=30, 
	selection=NULL, limit=-1, showClusterInfo=TRUE) {

    #check if table exists in db
    if (length(which(table==listGenDB(db))) == 0)
	stop("Could not find table in database")

    #check if table is of type NSV	
    meta<-dbReadTable(db$db,"metaData") #meta = table data in memory
    index<-which(meta$name==table)  #find index of table
    if (meta$type[index]!="NSV")
	stop("Not an NSV table")
    #get metadata about the table
    meta<-as.character(subset(metaGenDB(db),name==table)["annotation"])
    x<-unlist(strsplit(meta,";"))	
    window <-sub("window=","",x[3])
    overlap <- sub("overlap=","",x[4])
    word <-sub("word=","",x[5])
    last_window <-sub("last_window=","",x[6])	

    emm <- EMM(measure=measure,threshold=threshold)

    nSequences <- nSequences(db, rank, name)
    hierarchy <- getHierarchy(db, rank, name)

    if (limit != -1) nSequences <- min(nSequences,limit)

    #clusterinfo stores the last_clustering details
    clusterInfo<-list(nSequences)
    ### Kullback can not handle 0 counts!
    if(measure=="Kullback") d <- d + 1

    cat("GenModel: Creating model for ",rank,": ",name,"\n",sep="")

    if (is.null(selection)){
	i<-0
	total<-0
	while(i<nSequences){
	    #get 100 sequences at a time
	    d <-getSequences(db,rank,name,table,limit=c(i,100))
	    n <- length(d)
	    n <- min(n,100)
	    ids <- attr(d,"id")[1:n]
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
    } else if (!is.null(selection)){
	d<-getSequences(db, rank, name, table, limit=limit)
	ids <- attr(d,"id")[selection]
	d<-d[selection]
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
	    hierarchy=hierarchy, measure=measure, 
	    window=window, word=word, overlap=overlap, 
	    last_window=last_window)

    if (showClusterInfo) genModel$clusterInfo <- clusterInfo

    class(genModel) <- "GenModel"	
    genModel		
}

#	Returns the model states and the ID and the segment number of the sequences that are part of that state in format id:segment.
#	By default, returns all states as a list, if a modelState is specified. returns only the sequences that are part of that state
getClusteringDetails <- function(model, state=NULL)
{	
    
### FIXME: this does not work if a sequence has the same state twice!

    if(!is.null(state)) {
	occ <- lapply(model$clusterInfo, FUN=function(x) which(x==state))
	occ[sapply(occ, length) <1] <- NULL     
	return(data.frame(sequence=names(occ), segment=unlist(occ), 
			row.names=NULL))
    }
    
    l <- lapply(clusters(model$model), FUN=function(x) 
	    getClusteringDetails(model, state=x))

    names(l) <- clusters(model$model)
    l
}


#	Returns the sequences that are part of a given model state as list of DNA or NSV sequence objects
#
getClusteringSequences <- function(db, model, state, table="sequences")
{	
    #get the ids that are part of the model state as a list
    ids<-getClusteringDetails(model, state)

    #create different lists based on whether sequence or NSV table
    if (table =="sequences")
		stateSequences <- DNAStringSet()
    else
		stateSequences <- list()
    #loop through all the ids that are part of the modelstate		
    for(i in 1:nrow(ids))
    {
		id <- ids[i,]
		#split the sequence on ":" and get sequenceid and segment
		#sequence <- unlist(strsplit(id,split=":"))[1]
		sequence <- as.numeric(id["sequence"])
		#segment <- as.numeric(unlist(strsplit(id,split=":"))[2])
		segment <- as.numeric(id["segment"])
		#get window size from the model
		window <- as.numeric(model$window)
		#get the start position of the segment eg:1, 101, 201 etc
		start <- (segment - 1) * window + 1
		#get sequences and append to the list
		stateSequences <- c(stateSequences,getSequences(db, rank="id",name=sequence, table, start = start, length=window))
    }
    if(table!="sequences")
		class(stateSequences)<-"NSVSet"

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
    cat(x$rank,": ");cat(x$name,sep=", ");cat("\n")
    cat("\nModel:\n")
    print(x$model)
}

### plot a model
plot.GenModel <- function(x, ...) {
    plot(x$model,main=paste(x$rank,": ", x$name, sep=""), ...)
}



prune.GenModel <- function(x, ...) {
    prune(x$model, ...)

}

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
	limit=-1, ...) 
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
    rankNames <- getRank(db, rank)[,1]
    for(n in rankNames) {
	emm <- GenModelDB(db, table="NSV", rank, name=n,
		selection=selection, limit=limit, ...)
	n<-gsub("/","",n)
	saveRDS(emm, file=paste(rankDir, "/", n, ".rds", sep=''))
    }
}






