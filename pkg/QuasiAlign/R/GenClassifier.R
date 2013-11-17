
GenClassifier <- function(classifier_dir = NULL) {
    if(!.isGenClassifier(classifier_dir)) stop("You need to specify a directory with GenModels!")
    classifier_dir <- normalizePath(classifier_dir)
    
    structure(list(classifier_dir = classifier_dir), class="GenClassifier")
}

.isGenClassifier <- function(classifier_dir) 
  if(is.null(classifier_dir)) FALSE else file.exists(file.path(classifier_dir, "isGenClassifier"))

print.GenClassifier <- function(x, ...) {
    loc <- x$classifier_dir
    cat("GenClassifier\nLocation:" , loc, "\n") 
    cat("Rank models:",
	    paste(basename(list.dirs(loc))[-1], 
		    collapse=", "), "\n")
}

removeGenClassifier <- function(object) {
    ### first check if it looks like a RDP directory!
    if(!.isGenClassifier(object$classifier_dir)) stop("The given GenClassifier/directory does not look valid! Please remove the directory manually!")

    unlink(object$classifier_dir, recursive=TRUE)
}


### FIXME: add pruning option!

# Creates models in classifier_dir directory for all names in rank.  If selection is
# specified, then it uses only those sequences for creating the model
trainGenClassifier <- function(classifier_dir, db, rank = "Phylum", 
	table="NSV", selection=NULL, limit=NULL, add=FALSE, ...) 
{

  if(is(classifier_dir, "GenClassifier"))  
    classifier_dir  <- classifier_dir$classifier_dir
  else if(file.exists(classifier_dir) && !add) 
    stop("Directory/file exists! Set add=TRUE, remove old classifier, or choose a differnt directory!")
  
  ### check if classifier_dir exists
  if(!file.exists(classifier_dir)) {
    dir.create(classifier_dir)
    # mark as Genclassifier
    ff <- file(file.path(classifier_dir, "isGenClassifier"), "w")
    cat("TRUE", file=ff)
    close(ff)
  }
  
  ### create rank subdir
  rank <- tolower(rank)
  rankDir <- file.path(classifier_dir,rank)
  if(!file.exists(rankDir)) dir.create(rankDir)
  
  #get rank names
  if(is.null(selection))
    rankNames <- getRank(db, rank)
  else
    rankNames <- unique(rankNames_sel <- getRank(db, rank=rank , whereRank="id", 
                                                 whereName=selection))
  
  i <- 1  ### FIXME: otherwise R CMD check complains about missing global binding 
  foreach (i = 1:length(rankNames)) %dopar% {
    dbl <- reopenGenDB(db) ### reopen for multicore
    
    n <- rankNames[i]
    sel <- selection
    
    # only use selected sequences in model
    if(!is.null(selection)) {
      sel <- sel[rankNames_sel==n]  
      rn <- NULL
    } else rn <- n
    
    emm <- GenModelDB(dbl, rank=rank, name=rn, table=table,
                      selection=sel, limit=limit, ...)
    
    n <- gsub("/","",n)
    saveRDS(emm, file=paste(rankDir, "/", n, ".rds", sep=''))
    closeGenDB(dbl)
  }
  
  GenClassifier(classifier_dir)
}


# modelDir is a directory with subfolders for various ranks NSVList is a list
# containing NSV with a rank attribute and a "name" attribute which is a list
# of rankNames output is a data.frame containing the similarity scores,
# predicted value and the actual value
predict.GenClassifier <- function(object, newdata, rank="Phylum", 
	method="supported_transitions", match_cluster="exact", ...) {
#classify<-function(modelDir, NSVList, rank, method="supported_transitions")
#{

    modelDir <- object$classifier_dir
    NSVList <- newdata

    rankDir<-file.path(modelDir, tolower(rank))
    
    if (!file.exists(rankDir)) stop("Model directory ",rankDir," not found!")
    if (length(NSVList) ==0) stop("No sequence to classify against") 
    
    modelFiles <- dir(rankDir, full.names=TRUE)    
    modelNames <- sub(".rds", "", basename(modelFiles))

    i <- 1  ### FIXME: otherwise R check complains about missing global binding 
    classificationScores <- foreach (i =1:length(modelNames),.combine=cbind) %dopar% {
	cat("classify: Creating score matrix for", modelNames[i],"\n")
	model<-readRDS(modelFiles[i])
	sapply(NSVList, FUN = function(x) scoreSequence(model, 
			x, method=method, match_cluster=match_cluster, ...))
    }    
    colnames(classificationScores) <- modelNames
    
    winner <- apply(classificationScores, MARGIN=1, which.max)
    prediction <- factor(winner, levels=1:ncol(classificationScores), labels=colnames(classificationScores))

	list(scores=classificationScores, 
	    class=prediction)

}

