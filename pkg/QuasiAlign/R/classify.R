# Takes the sequences from a directory and splits them up into training and
# test sets.  Uses the training sequences to create models and stores them in
# the modelDir directory.  The pctTest is the fraction of sequences used for
# testing.
validateModels<-function(db, modelDir, rank="phylum", table="NSV", pctTest=0.1, method="supported_transitions", limit=NULL, numRanks=NULL, top=TRUE, measure="Manhattan", threshold=30, prune=TRUE, count_threshold=5)
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
    rankNames <- getRank(db, rank, count=TRUE, removeUnknown=TRUE)[,1]
	if (!is.null(numRanks))
	{
		if (top==TRUE)
			rankNames <- head(rankNames, numRanks)
    		else
			rankNames <- tail(rankNames, numRanks)
	}
	#testIds <- vector()
	#for (i in 1:length(rankNames))
    #clean up the RDP directories
	unlink("rdp/sequences/*")
	unlink("rdp/taxonomy/*")
	unlink("rdp/test/*")
	testIds<-foreach (i=1:length(rankNames), .combine=c) %dopar%
    {
		db_local <- reopenGenDB(db, flags=SQLITE_RO)
		#get number of sequences in the rank
		n <- nSequences(db_local,rank, rankNames[i])
		#how many sequences should we use
		if(is.null(limit)) limit <- n
		else limit <- min(n,limit)
		#create selection vector for this rank
		#train is the number of training cases
		train<-as.integer(pctTrain*limit)
		#test is the number of test cases
		test<-as.integer(pctTest*limit)
		#create an array of all sequences indices
		ids <- getRank(db_local,rank="id",whereRank=rank, whereName=rankNames[i])
		#get which indices are to be used for training
		#if (train <= 0) next;
		sampleIds <- sample(ids,limit)
		train <- sample(sampleIds,train)
		#remove these from the x indices
		sampleIds <- setdiff(sampleIds,train)
		#get which indices are to be used for testing
		test <- sample(sampleIds,test)
		if (length(train) > 10) {
			emm<-GenModelDB(db_local,measure=measure, threshold=threshold, table="NSV", rank, name=rankNames[i], selection=train)
			if (prune)
				emm <- prune(emm, count_threshold=count_threshold, transitions=TRUE)
		
		#save the model to file
		#some species names have "/" in them, need to remove them
		rankNames[i]<-gsub("/","",rankNames[i])
		saveRDS(emm, file=paste(rankDir, "/", rankNames[i], ".rds", sep=''))
		#takes the ids and creates a RDP training file for them
		createRDPTraining(db_local,rank,rankNames[i],train)
		closeGenDB(db_local)
		rm(db_local)
		if (length(test) > 0)
			test
		}
	}
 	
	testNames <- getRank(db, rank=rank, whereRank="id", whereName=testIds, all=TRUE, partialMatch=FALSE)
	testList <- getSequences(db, table="NSV", rank="id", name=testIds)
	#combine the rdp files
	rdpSequences <- list.files("rdp/sequences")
	for(i in 1:length(rdpSequences))
	{
		command <- paste("cat rdp/sequences/", rdpSequences[i]," >> rdp/sequences/train.fasta",sep="")
		system(command)
		unlink(paste("rdp/sequences/",rdpSequences[i],sep=""))
	}	
	rdpTaxonomy <- list.files("rdp/taxonomy")
	for(i in 1:length(rdpTaxonomy))
	{
		command <- paste("cat rdp/taxonomy/", rdpTaxonomy[i]," >> rdp/taxonomy/taxonomy.txt",sep="")
		system(command)
		unlink(paste("rdp/taxonomy/",rdpTaxonomy[i],sep=""))
	}

	testSequences <- getSequences(db, rank="id", name=testIds)
	#create test file in fasta format for RDP
	if (!file.exists("rdp/test"))
		dir.create("rdp/test", recursive=TRUE)
	write.XStringSet(testSequences, filepath="rdp/test/test.fasta")
	#run the RDP classifier
	#check if the macqiime or qiime exists
	if (Sys.which("macqiime") !="")
		exec <- "macqiime"
	else if (Sys.which("qiime") !="")
		exec <- "qiime"
	else
	{
		print("unable to find macqiime/qiime for RDP")
		exec <- NULL
	}
	if(!is.null(exec))
	{
		command <- "assign_taxonomy.py -i rdp/test/test.fasta -t rdp/taxonomy/taxonomy.txt -r rdp/sequences/train.fasta -o ."
		system(exec, input=command)
	}

	#	
	attr(testList,"rank")<-rank
    attr(testList,"name")<-testNames
	attr(testList,"id") <- testIds
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
    if (length(NSVList) ==0) stop("No sequence to classify against") 
    modelFiles <- dir(rankDir, full.names=TRUE)    
    modelNames <- sub(".rds", "", basename(modelFiles))
    #classificationScores <- matrix(NA, ncol=length(modelFiles), 
	#    nrow=length(NSVList))
    #colnames(classificationScores) <- modelNames

    classificationScores <- foreach (i =1:length(modelNames),.combine=cbind) %dopar% {
		cat("classify: Creating score matrix for", modelNames[i],"\n")
		model<-readRDS(modelFiles[i])
		#classificationScores[,i] <- sapply(NSVList, FUN =
		#	function(x) scoreSequence(model, x, method=method, plus_one=TRUE))
		sapply(NSVList, FUN =
			function(x) scoreSequence(model, x, method=method, plus_one=TRUE))
    }    
    
	colnames(classificationScores) <- modelNames
    winner <- apply(classificationScores, MARGIN=1, which.max)
	prediction <- matrix(NA,nrow=nrow(classificationScores), ncol=3)
    for(i in 1:nrow(classificationScores)) {
		prediction[i,1]<-rownames(classificationScores)[i]
		prediction[i,2]<-names(which.max(classificationScores[i,]))
		prediction[i,3] <- attr(NSVList,"name")[which(attr(NSVList,"id")==rownames(classificationScores)[i])]
	}
	colnames(prediction) <- c("id","predicted","actual")
	#prediction <- colnames(classificationScores)[winner]
	#prediction <- modelNames[winner]

    #actual <-  attr(NSVList, "name")
	#id <- attr(NSVList,"id")
    list(scores=classificationScores, 
	    prediction=prediction)
		#prediction=cbind(id=id,predicted=prediction, actual=actual))
}

createRDPTraining <- function(db, rank, name, ids)
{
	tax <- getHierarchy(db,rank="id",name=ids)
	#tax is a matrix
	outTrainTax <- vector()
	if (!file.exists("rdp/taxonomy"))
		dir.create("rdp/taxonomy", recursive=TRUE)
	if (!file.exists("rdp/sequences"))
		dir.create("rdp/sequences", recursive=TRUE)
	for(i in 1:nrow(tax))
	{
		header =""
		#seq = ""
		id = tax[,"Id"][i]
		seq = getSequences(db,rank="id",name=id,partialMatch=FALSE)
		header<- paste(id,"\t",sep="")
		header <- paste(header,"Root",sep="")
		for(j in 2:6)
		{
			header <- paste(header, ";" ,substr(tolower(names(tax[i,j])),1,1),"__",tax[i,j],sep="")
		}
		header <- gsub("\\(class\\)","",header)
		header <- gsub(" ","",header)		
		
		outTrainTax[i] <- header
		write.XStringSet(seq,filepath=paste("rdp/sequences/",rank,name,".fasta",sep=""), append=TRUE)
	}
	write(outTrainTax, file=paste("rdp/taxonomy/",rank,name,".txt",sep=""), append=TRUE)
}

