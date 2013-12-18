
### Constructor
NSVSet <- function(x=NULL, window=NULL, word=NULL, 
	last_window=NULL, overlap=NULL, rank=NULL, name=NULL) {
    structure(x, window=window, word=word, last_window=last_window, 
	    overlap=overlap, rank=rank, name=name, 
	    class="NSVSet")
}

### combine
c.NSVSet <- function(..., recursive = FALSE) {
    args <- list(...)
    window <- attr(args[[1]], "window")
    word <- attr(args[[1]], "word")

    args <- lapply(args, unclass)
    NSVSet(do.call("c", args), window=window, word=word)
}

### print an NSVSet object
print.NSVSet <- function(x, ...) {
    cat("Object of class NSVSet for", length(x), "sequences\n")
    cat("Window:", attr(x, "window"), "/ Word:", attr(x, "word"), "\n")
    #	    "/ last_window:", attr(x, "last_window"), "\n")
    #cat("Rank:", attr(x, "rank"), "/ Name:", attr(x, "name", exact=TRUE),"\n") 
    if(length(x)>0) {
	cat("Number of segments (table with counts):") 
	print(table(sapply(x, nrow)))
    }
}

### subset
`[.NSVSet` <- function(x, i, j, ..., drop = TRUE) { 
    r <- unclass(x)
    if(!missing(i)) r <- r[i] 
    if(!missing(j)) r <- lapply(r, "[", j, ,drop=FALSE) 

    NSVSet(r, window=attr(x, "window"), word=attr(x, "word"))
}

#plots a barplot with whiskers
plot.NSVSet <- function(x, ..., whiskers=TRUE)
{
    mean <- rowMeans(sapply(x, colMeans))
    minVal <- apply(sapply(x, apply, MARGIN=2, min), MARGIN=1, min)
    maxVal <- apply(sapply(x, apply, MARGIN=2, max), MARGIN=1, max)

    bp <- barplot(mean, ylim=c(0,max(maxVal)), las=2, ...)
    if(whiskers) errbar(bp, mean, maxVal, minVal, cap=0.005, add=T, pch=NA)
    invisible(bp)
}


### convert to NSVs
createNSVSet <- function(x, window=100, overlap=0, word=3, 
	last_window=FALSE, allOffsets=FALSE) {
    s <- lapply(Biostrings::as.list(x), FUN= function(i) .counter(i, window, 
		    overlap, word, last_window, allOffsets ))

    NSVSet(s, window=window, overlap=overlap, word=word, 
	    last_window=last_window)
}


createNSVTable <- function(db, tableFrom="sequences", table="NSV", rank=NULL, name=NULL, 
	window=100, overlap=0, word=3, last_window=FALSE, 
	limit=NULL, removeUnknownSpecies=FALSE, allOffsets=FALSE) {

    if(length(grep(" ", table))) stop("table cannot contain spaces!")
    
    if (any(table==listGenDB(db)))
	stop("A table with this name already exists in the db")
    
    NSV <- "id TEXT PRIMARY KEY REFERENCES classification(id), data BLOB"
    tr <- try(
	    dbSendQuery(db$db, 
		    statement = paste("CREATE TABLE ", table,
			    "(",  NSV,  ")", sep='')
		    )
	    )
    
    if(is(tr, "try-error")) stop("Unable to create table ", table)
    
    # insert into metadata
    meta<-paste("'", table, "','NSV','rank=", rank, ";name=",
	    name , ";window=", window, ";overlap=", overlap,
	    ";word=", word, ";last_window=", last_window , ";'", sep='')
    try(
	    dbSendQuery(db$db, 
		    statement = paste("INSERT INTO metaData VALUES ( ",  
			    meta,  ")", sep='')
		    )
	    )	
    
    ## number of sequences
    n <- nSequences(db, rank, name, table=tableFrom)
    if(n<1) {
	warning("No sequences to process!")
	return()
    }

    cat("CreateNSVTable: Number of sequences to process:", n, "\n")
    if(!is.null(limit) && limit<n) {
	cat("CreateNSVTable: Converting only", limit, "sequences (limit)\n")
	n <- limit
    } 

    ### process block sequences on each core
    block <- 100
    
    ### concurrent write does not work for SQLite
    ### so we write the result to file and then read it back
    tmp_dir <- paste(db$dbName, "_chunks", sep="")
    dir.create(tmp_dir)

    res <- foreach(start=seq(0, n-1, by=block), .combine=rbind) %dopar% {
#	ok <- 0
#	fail <- 0
#	total <- 0
	
	dbl <- reopenGenDB(db)
	
	s <- getSequences(dbl, rank, name, table=tableFrom,
		limit=c(start, block))
	n <- createNSVSet(s,  window, overlap, word, last_window, allOffsets)
	n <- lapply(n, FUN=function(x) base64encode(serialize(x, NULL)))

	n <- data.frame(id=I(names(s)), data=I(unlist(n)))
	saveRDS(n, file=file.path(tmp_dir, start))

	closeGenDB(dbl) ### this does not close the main db!
#	c(ok, fail, total)
	nrow(n)
    }
    #end loop
    
    ### write files sequential back to db 
    for(f in list.files(tmp_dir, full.names=TRUE)) {
	n <- readRDS(file=f)
	dbWriteTable(db$db, table, n, append = TRUE, row.names=FALSE)
    }
    unlink(tmp_dir, recursive=TRUE)
     
    res <- sum(res)

#    if(is.matrix(res)) res <- colSums(res)
#    names(res) <- c("ok", "fail", "total")
#    cat("CreateNSVTable: Read ", res["total"], " entries. Added ", 
#	    res["ok"] , " entries (",res["fail"]," failed).\n", sep="")

    cat("CreateNSVTable: Processed ", res, " sequences.\n", sep="") 
}


getNSVs <- function(db,  rank=NULL, name=NULL, 
	table="NSV", limit=NULL, random=FALSE, start=NULL, length=NULL, 
  removeUnknownSpecies=FALSE, annotation="id") {

    # limit = number of sequences to limit	
    # random = whether the sequences should be random
    # start = start of the chunk eg: 1
    # length = length of the chunk eg: 100 (should be called width?)
  
  if(!is.null(limit)) 
    limit <- paste(" LIMIT ", paste(limit,collapse=","))
  else limit <- "" 
  
  if(random) limit <- paste(" ORDER BY RANDOM() ", limit)
  
  if (!is.null(rank)) {    
    fullRank<-BioTools:::.pmatchRank(db, rank)
    #Do this so that the column order appears as 'order' since ORDER is a SQL keyword
    fullRankSQL<-paste("classification.'",fullRank,"'",sep="")
  }
  else fullRankSQL <-"-1"
  
  res <- dbGetQuery(db$db, statement = paste("SELECT ", table, 
                                             ".data AS data, classification.id AS id, ", 
                                             fullRankSQL ," AS fullRank  FROM sequences ", 
                                             " INNER JOIN classification ON classification.id = sequences.id INNER JOIN ", table, " ON ", table, ".id=sequences.id ",
                                             BioTools:::.getWhere(db, rank, name,
                                                                  removeUnknownSpecies), 
                                             limit, sep="")
  )

  if (nrow(res) == 0) stop("No rows found in the database")
  
  #get metadata about the table
  meta <- as.character(subset(metaGenDB(db),name==table)["annotation"])
  meta <- sapply(unlist(strsplit(meta,";")), strsplit, "=")	
  meta <- structure(sapply(meta, "[", 2), names=sapply(meta, "[", 1)) 
  
  ret <- lapply(res$data,decodeSequence)
  names(ret) <- res$id
  
  NSVSet(ret, 
         window=as.integer(meta["window"]),
         overlap=as.integer(meta["overlap"]), 
         word=as.integer(meta["word"]), 
         last_window=as.logical(meta["last_window"]), 
         rank=meta["rank"], 
         name=meta["name"])
}


### overwrite getSequences so it can decode NSVs
getSequences <- function(db,  rank=NULL, name=NULL,
                         table="sequences", limit=NULL, random=FALSE, start=NULL, 
                         length=NULL, removeUnknownSpecies=FALSE, 
                         annotation=Annotation_Id) {
  
  ### FIXME: check metadata table
  type <- metaGenDB(db, table)[,"type"]
  if(length(type) !=1) stop("table does not exist!")    
  
  if (type=="sequence") getX <- BioTools::getSequences
  else getX <- getNSVs
  
  getX(db,  rank, name,
       table, limit, random, start, 
       length, removeUnknownSpecies, 
       annotation)
}

# reads all fasta files in a directory into a db and 
# creates NSV table with all sequences
processSequences <- function(dir, db, annotation=Annotation_Greengenes,
	...) {
    for(f in dir(dir, full.names=TRUE))
    {
	cat("Processing file: ",f,"\n")
	addSequences(db, f, annotation=annotation)
    }

    createNSVTable(db, ...)
}

## helper
decodeSequence <- function(sequence) {
    if(length(sequence)==1) unserialize(base64decode(sequence, what="raw"))
    else lapply(sequence, FUN=function(x) unserialize(base64decode(x, what="raw")))
    }




