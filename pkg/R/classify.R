# Takes the sequences from a directory and splits them up into training and
# test sets.  Uses the training sequences to create models and stores them in
# the modelDir directory.  The pctTest is the fraction of sequences used for
# testing.
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
	rankNames[,1][i]<-gsub("/","",rankNames[,1][i])
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

# modelDir is a directory with subfolders for various ranks NSVList is a list
# containing NSV with a rank attribute and a "name" attribute which is a list
# of rankNames output is a data.frame containing the similarity scores,
# predicted value and the actual value
classify<-function(modelDir, NSVList)
{

    ### takes models from rank subdir
    rank <- tolower(attr(NSVList,"rank"))
    
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
	cat("classify: Creating score matrix for", modelName,"\n")
	modelName<- sub(".rds","",modelName)
	model<-readRDS(modelFile)
	modelSim<-data.frame()
	for (NSV in NSVList)
	{
	    sc<-score(model$model,NSV,plus_one=TRUE)
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


