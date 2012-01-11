
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

	rank <- .pmatchRank(db, rank)
	name<- unlist(attr(d,"name"))
	op <- paste(rank,": ", name)	
	seq <- paste(length(d)," sequences")
	op<- c(op, seq )
	
	genModel <- list(name=name, rank=rank, model=emm)
	class(genModel) <- "genModel"
	
	genModel	
	
}


processSequencesGreengenes <- function(dir, db) {
	for(f in dir(dir, full.names=T))
	{
	    cat("Processing file: ",f,"\n")
	    addSequencesGreengenes(db, f)
	}
	    
	createNSVTable(db, "NSV")
}


### validateModels...
    # create a selection information
    # call createModels
    # classify (with data from db selection information)
    

## sel <- sample(c(0,1), 1000, prob=c(.9,.1), replace=TRUE)
## sel <- sample(c(rep(1, times=100), rep(0, times=900)))



## FIXME: Test/Training selection
createModels <- function(modelDir, rank = "phylum", db) # selection argument
{
    ### check if modelDir exists
    ### create rank subdir
    
	rankNames <- getRank(db, rank)[,1]
	for(n in rankNames) {
	    emm <- genModel(db, table="NSV", rank, name=n)
	    cat("Creating model for ", rank, ":", n, "\n")
	    
	    saveRDS(emm, file=paste(modelDir, "/", n, ".rds", sep=''))
	}
}



## pass a list of NSVs instead of file.
## need a processSequenceFile (not database)

classify<-function(modelDir, seqFile, rank)
{

     ### takes models from rank subdir

	### read and NSVs for test Seq

	modelFiles <- dir(modelDir, full.names=TRUE)    
	    
	### for each model
	##	    score
	###	    cbind to bind to resulting data.frame


	### for each row in data.frame find max.
    ### output: return data.frame
    
	
	
	
	#outfile is where the results are saved - can make this a parameter
	outfile<-"ClassifyResults.txt"
	if (file.exists(outfile))
		unlink(outfile)
	

	modelNames<-vector()
	#header for the result file
	cat("Sequence\t",file=outfile,append=TRUE)
	for(f in dir(modelDir, full.names=F))
	{
	  f<-sub(".rds","",f)
	  modelNames<-c(modelNames,f)
	  cat(f,"\t",file=outfile,append=TRUE)
	}
	cat("Actual\tPredicted",file=outfile,append=TRUE)
	cat("\n",file=outfile,append=TRUE)
	#delete db file if already exists
	if (file.exists(".classify.sqlite"))
		unlink(".classify.sqlite")
	#create new temp db for storing seqFile
	.dbc<-createGenDB(".classify.sqlite")
	#read file
	addSequencesGreengenes(.dbc,seqFile)
	#convert sequences to NSV
	createNSVTable(.dbc,"NSV")
	#get all sequences as a list
	sequences<-getSequences(.dbc,table="NSV")
	#get all the values in the rank
	rankValues<-getHierarchy(.dbc,rank)[,1]
	#loop through all the sequences
	for(i  in 1:length(sequences))
	{
	 #loop through each model in the modelDir and find similarity	 
		cat(i,"\t",file=outfile,append=TRUE)
		maxPosition<-0
		maxValue <- 0
		j<-1
		for(f in dir(modelDir, full.names=T))
		{
		  .model<-readRDS(f)
		  score<-score(.model$model,sequences[[i]]+1,plus_one=TRUE)
		  if (score > maxValue) {
			maxPosition<-j
			maxValue<-score
		  }
		  j<- j+1
		  cat(score,"\t",file=outfile,append=TRUE)
		}
		#output the actual rank
		cat(rankValues[i],"\t",file=outfile,append=TRUE)
		#output predicted rank
		cat(modelNames[maxPosition],"\t",file=outfile,append=TRUE)
		cat("\n",file=outfile,append=TRUE)
	}
	unlink(".classify.sqlite")
	cat("Done. Results are in file ",outfile," \n")
	
}
