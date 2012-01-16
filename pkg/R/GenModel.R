
genModel <- function(db, rank=NULL, name=NULL, table, limit=-1, 
	measure="Kullback", threshold=0.10, plus_one=TRUE, selection=NULL) {
	
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
    
validateModels<-function(dir, modelDir, rank="phylum", pctTrain=0.9, pctTest=0.1, db)
{
	#dir => directory containing FASTA files which are to be used for model
	#modelDir => directory where models are to be stored
	#pctTrain => percentage of each rank to be used for creating the training model
	#pctTest => percentage of each rank to be used for testing the model
	
	
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
	#read sequences and convert to NSV
	processSequencesGreengenes(dir, db)
	#create a list with a vector of selection for EACH rank
	trainingList<-list()
	testList<-list()
	testNames<-list()
	#get all the  rankNames for the given rank
	rankNames <- getRank(db, rank)

	for (i in 1:length(rankNames[,1]))
	{
		#get number of sequences in the rank
		n<- nSequences(db,rank, rankNames[,1][i])
		#create selection vector for this rank
		#train
		train<-as.integer(pctTrain*n)
		test<-as.integer(pctTest*n)
		#create a list with selected=0 for all initially
		selList<- list(val=as.vector(1:n),training=as.vector(rep(0,n)))
		#samp contains all the sequences which have been selected for training
		samp<- sample(selList$val,size=train,replace=FALSE)
		#mark training as true for all the samp
		selList$training[samp]<-1
		#get indices which have been selected for training
		sel<-selList$val[selList$training==1]
		#get the remaining indices which are for testing
		notsel<-selList$val[selList$training==0]
		#append selection vector to main list
		trainingList<-c(trainingList, sel)
		#create model using the training set
		emm<-genModel(db, table="NSV", rank, name=rankNames[,1][i], selection=sel)
	    #save the model to file
		saveRDS(emm, file=paste(rankDir, "/", rankNames[,1][i], ".rds", sep=''))
		#########
		# can also do this with the createModels function as:
		# createModels(modelDir, rank, db, sel) => but need to find way to add rankName
		##########
		#get all sequences and filter it to just test sequences
		d<-getSequences(db,table="NSV",rank=rank,name=rankNames[,1][i])
		#filter sequences and add to test list
		testList<-c(testList,d[notsel])
		#Names are lost after filtering, so need to keep a list of ranknames
		testNames<-c(testNames,attr(d,"name")[[1]][notsel])
	}
	#add attributes to test
	attr(testList,"rank")<-rank
	attr(testList,"name")<-testNames
	classify(modelDir,testList)
}
## sel <- sample(c(0,1), 1000, prob=c(.9,.1), replace=TRUE)
## sel <- sample(c(rep(1, times=100), rep(0, times=900)))

# Creates models in modelDir directory for all names in rank.
# If selection is specified, then it uses only those sequences for creating the model
createModels <- function(modelDir, rank = "phylum", db, selection=NULL) 
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
	    emm <- genModel(db, table="NSV", rank, name=n,selection=selection)
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
			actualValues<-rbind(attr(NSVList,"name")[[i]], deparse.level=3)
		else 
			actualValues<-rbind(actualValues,attr(NSVList,"name")[[i]], deparse.level=3)
	}
	colnames(predValues)<-"predicted"
	colnames(actualValues)<-"actual"
	classificationScores<-cbind(classificationScores,actualValues)
	classificationScores<-cbind(classificationScores,predValues)
	return(classificationScores)
	
}
