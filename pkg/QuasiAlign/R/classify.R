# Takes the sequences from a directory and splits them up into training and
# test sets.  Uses the training sequences to create models and stores them in
# the modelDir directory.  The pctTest is the fraction of sequences used for
# testing.
validateModels<-function(db, modelDir, rank="phylum", table="NSV", pctTest=0.1, method="supported_transitions", limit=NULL, numRanks=NULL, top=TRUE, measure="Manhattan", threshold=30, count_threshold=5)
{
    #dir => directory containing FASTA files which are to be used for model
    #modelDir => directory where models are to be stored
    #pctTrain => percentage of each rank to be used for creating the training model
    #pctTest => percentage of each rank to be used for testing the model
    #db=> database where sequences are to be stored	
    #numRanks=> number of names of each rank to consider , for example just the top 5 phylums, etc
    rankDir<-file.path(modelDir,rank)
    if (file.exists(modelDir))
    {
		if(!file.exists(rankDir)) dir.create(rankDir)
		else
			#remove all existing model files
			unlink(file.path(rankDir,dir(rankDir)))
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
    rankNames <- getRank(db, rank, count=TRUE, removeUnknown=TRUE)
	if (!is.null(numRanks))
	{
		if (top==TRUE)
			rankNames <- head(rankNames, numRanks)
    		else
			rankNames <- tail(rankNames, numRanks)
	}
    testIds<-foreach (i=1:length(rankNames[,1])) %dopar%
    {
		db_local <- reopenGenDB(db, flags=SQLITE_RO)
		#get number of sequences in the rank
		n <- nSequences(db_local,rank, rankNames[,1][i])
		#how many sequences should we use
		if(is.null(limit)) limit <- n
		else limit <- min(n,limit)
		#create selection vector for this rank
		#train is the number of training cases
		train<-as.integer(pctTrain*limit)
		#test is the number of test cases
		test<-as.integer(pctTest*limit)
		#create an array of all sequences indices
		ids <- getRank(db_local,rank="id",whereRank=rank, whereName=rankNames[,1][i])[,1]
		#get which indices are to be used for training
		#if (train <= 0) next;
		sampleIds <- sample(ids,limit)
		train <- sample(sampleIds,train)
		#remove these from the x indices
		sampleIds <- setdiff(sampleIds,train)
		#get which indices are to be used for testing
		test <- sample(sampleIds,test)
		if (length(train) > 0) {
			emm<-GenModelDB(db_local,measure=measure, threshold=threshold, table="NSV", rank, name=rankNames[,1][i], selection=train)
			emm <- prune(emm, count_threshold=count_threshold)
		}
		#save the model to file
		#some species names have "/" in them, need to remove them
		rankNames[,1][i]<-gsub("/","",rankNames[,1][i])
		saveRDS(emm, file=paste(rankDir, "/", rankNames[,1][i], ".rds", sep=''))
		closeGenDB(db_local)
		rm(db_local)
		if (length(test) > 0)
			test
		
		#get all sequences and filter it to just test sequences
		#d<-getSequences(db,table="NSV",rank=rank,name=rankNames[,1][i])
		#filter sequences and add to test list
		#testList <- c(testList,d[test])
		#testNames <- c(testNames, attr(d,"name")[test])	
    }
    #if(length(testList)==0)
	#	stop("Insufficient sequences have been selected for testing")
    #add attributes to test
	testIds <- unlist(testIds)
	testNames <- getRank(db, rank=rank, whereRank="id", whereName=testIds, partialMatch=FALSE)[,1]
	testList <- getSequences(db, table="NSV", rank="id", name=testIds)
	attr(testList,"rank")<-rank
    attr(testList,"name")<-testNames
	if (length(method) == 1)
    	return(classify(modelDir, testList, rank=rank, method=method))
	else if (length(method) > 1)
	{
		ret<-list()
		for(i in 1:length(method))
				ret[[i]]<-classify(modelDir, testList, rank=rank, method=method[i])
		return(ret)
	}
}


# modelDir is a directory with subfolders for various ranks NSVList is a list
# containing NSV with a rank attribute and a "name" attribute which is a list
# of rankNames output is a data.frame containing the similarity scores,
# predicted value and the actual value
classify<-function(modelDir, NSVList, rank, method="supported_transitions")
{

    rankDir<-file.path(modelDir, tolower(rank))
    if (!file.exists(rankDir)) stop("Model directory ",rankDir," not found!")
    
    modelFiles <- dir(rankDir, full.names=TRUE)    
    modelNames <- sub(".rds", "", basename(modelFiles))
    #classificationScores <- matrix(NA, ncol=length(modelFiles), 
	#    nrow=length(NSVList))
    #colnames(classificationScores) <- modelNames

    classificationScores <- foreach (i =1:length(modelFiles),.combine=cbind) %dopar% {
	
	cat("classify: Creating score matrix for", modelNames[i],"\n")
	model<-readRDS(modelFiles[i])
	#classificationScores[,i] <- sapply(NSVList, FUN =
	#	function(x) scoreSequence(model, x, method=method, plus_one=TRUE))
	sapply(NSVList, FUN =
		function(x) scoreSequence(model, x, method=method, plus_one=TRUE))
    }    
    
	colnames(classificationScores) <- modelNames
    winner <- apply(classificationScores, MARGIN=1, which.max)
    prediction <- modelNames[winner]

    actual <-  attr(NSVList, "name")

    list(scores=classificationScores, 
	    prediction=cbind(predicted=prediction, actual=actual))
}


