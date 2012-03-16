# creates an model from sequences in the db 
genModel <- function(db, rank=NULL, name=NULL, table, 
	measure="Manhattan", threshold=30, plus_one=TRUE, 
	selection=NULL, limit=-1) {

    #check if table exists in db
    if (length(which(table==listGenDB(db))) == 0)
		stop("Could not find table in database")
    #check if table is of type NSV	
    meta<-dbReadTable(db$db,"metaData") #meta = table data in memory
    index<-which(meta$name==table)  #find index of table
    if (meta$type[index]!="NSV")
		stop("Not an NSV table")
    emm <- EMM(measure=measure,threshold=threshold)
    # how many sequences
	nSequences<-nSequences(db,rank,name)
    # if limit has been specified, use that
	if (limit != -1)
		nSequences <- min(nSequences,limit)
    cat("genModel: Creating model for ",rank,": ",name,"\n",sep="")
    if (is.null(selection))
    {
		i<-0
		total<-0
		while(i<nSequences)
		{
	    	d<-getSequences(db,rank,name,table,limit=c(i,100))
	    	i<-min(i+100,nSequences)
	    	sequences<-.make_stream(d)

	    	### Kullback can not handle 0 counts!
	    	if(measure=="Kullback") sequences<-sequences + 1

	    	build(emm,.make_stream(d)+1)
	    	reset(emm)
	    	cat("genModel: Processed",i,"sequences\n")
		}
    }
    else if (!is.null(selection))
    {
		d<-getSequences(db, rank, name, table, limit=limit)
		d<-d[selection]
		if (length(d)==0)
	    	stop("GenModel called with 0 sequences")
		for(i in 1:length(d))
		{
	    	sequence<- d[[i]]		

	    	if(measure=="Kullback") sequence <- sequence + 1

	    	build(emm,sequence)
	    	reset(emm)
	    	if (i%%100==0) cat("genModel: Processed",i,"sequences\n")
		}	
		cat("genModel: Processed",length(d),"sequences\n")
    }
    rank <- .pmatchRank(db, rank)
    rankName<- unique(unlist(attr(d,"name")))
    nSequences<-length(d)	
    genModel <- list(name=rankName, rank=rank, nSequences=nSequences, model=emm)
    class(genModel) <- "GenModel"	
    genModel		
}

### print basic info about a model
print.GenModel <- function(x, ...) {
    cat("Object of class GenModel with", x$nSequences, "sequences\n")
    cat(x$rank,": ", x$name, "\n", sep="")

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
processSequences <- function(dir, db, reader = addSequencesGreengenes, 
	...) {
    for(f in dir(dir, full.names=TRUE))
    {
	cat("Processing file: ",f,"\n")
	reader(db, f)
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
	emm <- genModel(db, table="NSV", rank, name=n,
		selection=selection, limit=limit, ...)
	n<-gsub("/","",n)
	saveRDS(emm, file=paste(rankDir, "/", n, ".rds", sep=''))
    }
}






