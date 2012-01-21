
genModel <- function(db, rank=NULL, name=NULL, table, 
	measure="Kullback", threshold=0.10, plus_one=TRUE, selection=NULL, limit=-1) {
	
	### FIXME: check in metadata if it is NSV

	emm <- EMM(measure=measure,threshold=threshold)
	d<-getSequences(db, rank, name, table, limit=limit)
	if (!is.null(selection))
		d<-d[selection]
	else if (length(d)==0)
		stop("GenModel called with 0 sequences")
	
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

# reads all fasta files in a directory into a db and 
# creates NSV table with all sequences
processSequencesGreengenes <- function(dir, db) {
	for(f in dir(dir, full.names=TRUE))
	{
	    cat("Processing file: ",f,"\n")
	    addSequencesGreengenes(db, f)
	}
	    
	createNSVTable(db, "NSV")
}

sequencesToModels <- function(dir, modelDir, rank) {
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
	for(f in dir(dir, full.names=TRUE))
	{
	    cat("Processing file: ",f,"\n")
		fileName <-basename(f)
		modelName<-sub(".fasta",".rds",fileName)
		dbName<-sub(".fasta",".sqlite",fileName)
		db<-createGenDB(dbName)
	    addSequencesGreengenes(db, f)
		createNSVTable(db, "NSV")
	    emm <- genModel(db, table="NSV")
	    saveRDS(emm, file=file.path(rankDir,modelName))
		unlink(dbName)
		rm(db)
		rm(emm)
	}
	    
}

# Takes the sequences from a directory and splits them up into training and test sets.
# Uses the training sequences to create models and stores them in the modelDir directory.
# The pctTest is the fraction of sequences used for testing.
validateModels<-function(db, modelDir, rank="phylum", table="NSV", pctTest=0.1)
{
	#dir => directory containing FASTA files which are to be used for model
	#modelDir => directory where models are to be stored
	#pctTrain => percentage of each rank to be used for creating the training model
	#pctTest => percentage of each rank to be used for testing the model
    #db=> database where sequences are to be stored	
	
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
	#db<-createGenDB(".validateModels.sqlite")
	#read sequences and convert to NSV
	#processSequencesGreengenes(dir, db)
	pctTrain = 1 - pctTest
	#create a list with a vector of selection for EACH rank
	trainingList<-list()
	testList<-list()
	testNames<-vector()
	#get all the  rankNames for the given rank
	rankNames <- getRank(db, rank)

	for (i in 1:length(rankNames[,1]))
	{
		#get number of sequences in the rank
		n <- nSequences(db,rank, rankNames[,1][i])
		#create selection vector for this rank
		#train
		train<-as.integer(pctTrain*n)
		test<-as.integer(pctTest*n)
		selList<-sample(c(rep(1,train),rep(0,test)))
		sel<-which(selList==1)
		notsel<-which(selList==0)
		#create model using the training set
		emm<-genModel(db, table="NSV", rank, name=rankNames[,1][i], selection=sel)
	    #save the model to file
		saveRDS(emm, file=paste(rankDir, "/", rankNames[,1][i], ".rds", sep=''))
		#get all sequences and filter it to just test sequences
		d<-getSequences(db,table="NSV",rank=rank,name=rankNames[,1][i])
		#filter sequences and add to test list
		testList<-c(testList,d[notsel])
		#Names are lost after filtering, so need to keep a list of ranknames
		#testNames[i,1]<-attr(d,"name")[[1]][notsel]
		testNames<-c(testNames,attr(d,"name")[notsel])
	}
	#add attributes to test
	unlink(".validateModels.sqlite")
	rm(db)
	#testNamesdf<-data.frame()
	#for(j in 1:length(testNames))
	#	testNamesdf[j,1]<-testNames[j]
	attr(testList,"rank")<-rank
	attr(testList,"name")<-testNames
	return(classify(modelDir,testList))
}
## sel <- sample(c(0,1), 1000, prob=c(.9,.1), replace=TRUE)
## sel <- sample(c(rep(1, times=100), rep(0, times=900)))

# Creates models in modelDir directory for all names in rank.
# If selection is specified, then it uses only those sequences for creating the model
createModels <- function(modelDir, rank = "phylum", db, selection=NULL, limit=-1) 
{
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
	    emm <- genModel(db, table="NSV", rank, name=n,selection=selection,limit=limit)
	    cat("Creating model for ", rank, ":", n, "\n")	    
	    saveRDS(emm, file=paste(rankDir, "/", n, ".rds", sep=''))
	}
}



# modelDir is a directory with subfolders for various ranks
# NSVList is a list containing NSV with a rank attribute and a "name" attribute which is a list of rankNames
# output is a data.frame containing the similarity scores, predicted value and the actual value
classify<-function(modelDir, NSVList)
{

     ### takes models from rank subdir
	rank<-attr(NSVList,"rank")
	if (is.null(rank))
		stop("NSVList does not have a rank attribute")
    rankDir<-file.path(modelDir,rank)
	### read and NSVs for test Seq
	if (!file.exists(rankDir))
		stop("Directory ",rankDir," not found")
	modelFiles <- dir(rankDir, full.names=TRUE)    
	classificationScores<-data.frame()
	for (modelFile in modelFiles)
	{
		modelName<-basename(modelFile)
		modelName<- sub(".rds","",modelName)
		model<-readRDS(modelFile)
		modelSim<-data.frame()
		for (NSV in NSVList)
		{
			sc<-score(model$model,NSV+1,plus_one=TRUE)
			if(length(modelSim)==0)
				modelSim<-rbind(sc,deparse.level=3) #deparse.level=3=>no rownames
			else
				modelSim<-rbind(modelSim,sc,deparse.level=3) #deparse.level=3=> no rownames
		}
		colnames(modelSim)<-modelName
		if (length(classificationScores)!=0)
			classificationScores<-cbind(classificationScores,modelSim)
		else if (length(classificationScores)==0)
			classificationScores <- modelSim
	}    
	#find predicted value
	predValues<-data.frame()
	actualValues<-data.frame()
	for (i in 1:length(classificationScores[,1]))
	{
		#find the max row
		maxRow=which.max(classificationScores[i,])
		predicted=colnames(classificationScores)[maxRow]
		if(length(predValues)==0)
			predValues<-rbind(predicted, deparse.level=3)
		else
			predValues<-rbind(predValues,predicted, deparse.level=3)
		if(length(actualValues)==0)
			#actualValues<-rbind(attr(NSVList,"name")[i,1], deparse.level=3)
			actualValues<-rbind(attr(NSVList,"name")[[i]], deparse.level=3)
		else 
			#actualValues<-rbind(actualValues,attr(NSVList,"name")[i,1], deparse.level=3)
			actualValues<-rbind(actualValues,attr(NSVList,"name")[[i]], deparse.level=3)
	}
	colnames(predValues)<-"predicted"
	colnames(actualValues)<-"actual"
	classificationScores<-cbind(classificationScores,actualValues)
	classificationScores<-cbind(classificationScores,predValues)
	return(classificationScores)
	
}
