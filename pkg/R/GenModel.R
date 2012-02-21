# creates an model from sequences in the db 
genModel <- function(db, rank=NULL, name=NULL, table, 
	measure="Kullback", threshold=0.10, plus_one=TRUE, 
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
	nSequences<-nSequences(db,rank,name)
	if (limit != -1)
		nSequences <- min(nSequences,limit)
	#d<-getSequences(db, rank, name, table, limit=limit)
	#if (!is.null(selection))
	#	d<-d[selection]
	#else if (length(d)==0)
	#	stop("GenModel called with 0 sequences")
    cat("genModel: Creating model for rank:",rank,",name:",name,"\n")
	if (is.null(selection))
	{
		i<-0
		total<-0
		while(i<nSequences)
		{
			d<-getSequences(db,rank,name,table,limit=c(i,100))
			i<-min(i+100,nSequences)
			cat("genModel: Read ",i,"sequences \n")
			sequences<-.make_stream(d)
			if(plus_one) sequences<-sequences + 1
			build(emm,.make_stream(d)+1)
			reset(emm)
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
			if (i%%100==0) cat("genModel: Read ",i," sequences \n")
			sequence<- d[[i]]		
			if(plus_one) sequence <- sequence + 1
			build(emm,sequence)
			reset(emm)
		}	
		cat("genModel: Read ",length(d)," sequences \n")
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
    cat("Object of class GenModel")
    cat(" of", x$nSequences, "sequences.\n")
    cat("Rank/Name:", x$rank,"-", x$name, "\n")
    
    cat("\nModel:\n")
    print(x$model)
}

### plot a model
plot.GenModel <- function(x, ...) {
    plot(x$model,main=paste(x$rank,"-", x$name), ...)
}

### score a new sequence of NSVs against a model
score.GenModel <- function(x, newdata, method = "prod", 
	match_cluster = "nn", plus_one=TRUE, 
	initial_transition = FALSE) {
    score(x$model, newdata=newdata, method=method, 
	    match_cluster=match_cluster, plus_one=plus_one, 
	    initial_transition=initial_transition)
}

# reads all fasta files in a directory into a db and 
# creates NSV table with all sequences
processSequences <- function(dir, db, reader = addSequencesGreengenes) {
	for(f in dir(dir, full.names=TRUE))
	{
	    cat("Processing file: ",f,"\n")
	    reader(db, f)
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
		rankName<-sub(".fasta","",fileName)
		#donot overwrite existing model files
		if (!file.exists(file.path(rankDir,modelName)))
		{
			if (file.exists(dbName))
				unlink(dbName)	
			db<-createGenDB(dbName)
    		addSequencesGreengenes(db, f)
			createNSVTable(db, "NSV")
    		emm <- genModel(db, table="NSV",rank=rank,name=rankName)
    		saveRDS(emm, file=file.path(rankDir,modelName))
			unlink(dbName)
			rm(db)
			rm(emm)
		}
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
		#sel contains the sequences which have been selected for training
		sel<-which(selList==1)
		#notsel contains the test sequences
		notsel<-which(selList==0)
		#create model using the training set
		emm<-genModel(db, table="NSV", rank, name=rankNames[,1][i], selection=sel)
	    #save the model to file
		saveRDS(emm, file=paste(rankDir, "/", rankNames[,1][i], ".rds", sep=''))
		#get all sequences and filter it to just test sequences
		d<-getSequences(db,table="NSV",rank=rank,name=rankNames[,1][i])
		#filter sequences and add to test list
		testList<-c(testList,d[notsel])
		testNames<-c(testNames,attr(d,"name")[notsel])
	}
	if(length(testList)==0)
		stop("Insufficient sequences have been selected for testing")
	#add attributes to test
	rm(db)
	attr(testList,"rank")<-rank
	attr(testList,"name")<-testNames
	return(classify(modelDir,testList))
}

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
		cat("clasify: Creating score matrix for: ",modelName,"\n")
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
	prediction<-data.frame() #final prediction
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
	prediction<-cbind(actualValues)
	prediction<-cbind(prediction,predValues)
	classification<-list(scores=classificationScores,prediction=prediction)
	return(classification)
	
}
