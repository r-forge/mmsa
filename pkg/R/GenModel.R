
genModel <- function(db, rank=NULL, name=NULL, table, limit=-1, 
	measure="Kullback", threshold=0.10, plus_one=TRUE) {
	
	### FIXME: check in metadata if it is NSV

	emm <- EMM(measure=measure,threshold=threshold)	
	d<-getSequences(db, rank, name, table, limit=limit)
	
	for(i in 1:length(d))
	{
		sequence<- d[[i]]		
		if(plus_one) sequence <- sequence + 1
		build(emm,sequence)
		reset(emm)
	}	

	# FIXME: add name attribute. getSequences should pass on attributes
	## for name and rank...
	rank <- .pmatchRank(db, rank)
	name<- unlist(attr(d,"name"))
	#op <- paste(rank,": ", name, " - ", length(d), " sequences" , sep = '')
	op <- paste(rank,": ", name)	
	op<- c(op, length(d), " sequences" )
	genModel <- list(name=op, model=emm)
	class(genModel) <- "genModel"
	
	genModel	
	
}

plot.genModel<-function(x, ...)
{
	plot(x$model, main=x$name,...)
	
}



processSequencesGreengenes <- function(dir, db) {
	for(f in dir(dir, full.names=T))
	{
	    addSequencesGreengenes(db, f)
	}
	    
	createNSVTable(db, "NSV")
}

createModels <- function(models, rank = "phylum", db)
{
	rankNames <- getRank(db, rank)[,1]
	for(n in rankNames) {
	    emm <- genModel(db, table="NSV", rank="phylum", name=n)
	    cat("Creating model for ", rank, ":", n, "\n")
	    
	    saveRDS(emm, file=paste(models, "/", n, ".rds", sep=''))
	}
}

classify<-function(sequence, models)
{
}
