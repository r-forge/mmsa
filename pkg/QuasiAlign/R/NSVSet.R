
### Constructor
NSVSet <- function(x=NULL, window=NULL, word=NULL, 
	last_window=NULL, overlap=NULL, rank=NULL, name=NULL) {
    structure(x, window=window, word=word, last_window=last_window, 
	    overlap=overlap, rank=rank, name=name ,class="NSVSet")
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
    ### This should work but as.list does not inside the package!
    #s <- lapply(x, .counter, window, overlap, word, last_window)
    s <- lapply(1:length(x), FUN= function(i) .counter(x[[i]], window, 
		    overlap, word, last_window, allOffsets ))

    NSVSet(s, window=window, overlap=overlap, word=word, 
	    last_window=last_window)
}


createNSVTable <- function(db, table="NSV", 
	rank=NULL, name=NULL, window=100,
	overlap=0, word=3, last_window=FALSE, limit=NULL, removeUnknownSpecies=FALSE, allOffsets=FALSE) {

    if(length(grep(" ", table))) stop("table cannot contain spaces!")
    if (length(which(table==listGenDB(db))) > 0)
	stop("A table with this name already exists in the db")
    NSV <- "id TEXT PRIMARY KEY REFERENCES classification(id), data BLOB"
    try(
	    dbSendQuery(db$db, 
		    statement = paste("CREATE TABLE ", table,
			    "(",  NSV,  ")", sep='')
		    )
	    )
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
    ok <- 0
    fail <- 0
    total <- 0

    #start loop
    start<-0 #start position of query
    num_records<-min(limit,100) # number of records at a time
    if (is.null(limit))
	limit <- 500000
    while(TRUE && total <limit)
    {
	if (is.null(rank) && is.null(name))
	{

	    d <- dbGetQuery(db$db, statement = 
		    paste("SELECT s.* FROM sequences s ", if (removeUnknownSpecies) " INNER JOIN classification c ON c.id=s.id WHERE c.[Species] NOT LIKE 'Unknown%' " ,"LIMIT ", start,",
			    ",num_records,sep=""))    
			}
			else {

			    #statement<-	paste("SELECT sequences.id, sequences.data FROM sequences INNER JOIN classification ON sequences.id=classification.id WHERE classification.",
			    #	.pmatchRank(db,rank)," LIKE '", 
			    #	name,"%'", if(removeUnknownSpecies) " AND classification.[Species] NOT LIKE 'Unknown%' ", " LIMIT ", start,", ",num_records, sep="")
			    statement<-	paste("SELECT sequences.id, sequences.data FROM sequences INNER JOIN classification ON sequences.id=classification.id ", BioTools:::.getWhere(db, rank, name) ,
				    if(removeUnknownSpecies) " AND classification.[Species] NOT LIKE 'Unknown%' ", " LIMIT ", start,", ",num_records, sep="")


				d <- dbGetQuery(db$db, statement = 
				    paste("SELECT sequences.id, sequences.data FROM sequences INNER JOIN classification ON sequences.id=classification.id ", BioTools:::.getWhere(db, rank, name) ,
					    if(removeUnknownSpecies) " AND classification.[Species] NOT LIKE 'Unknown%' ", " LIMIT ", start,", ",num_records, sep=""))    
			    }
			if (nrow(d)==0) break;

			dbBeginTransaction(db$db)
			for(i in 1:nrow(d))  {

			    #make NSV
			    nsv <- .counter(d$data[i], window, overlap, word, last_window, allOffsets)
			    d$data[i] <- base64encode(serialize(nsv, NULL))
			    #end make NSV
			    ## Insert into DB
			    dat <- paste("'",d[i,],"'", sep='', collapse=', ')

			    tr<- try(
				    dbSendQuery(db$db,          
					    statement = paste("INSERT INTO ", table,
						    " VALUES (", dat, ")", sep=''))
				    )

			    if(!is(tr, "try-error")) ok <- ok+1
			    else 
			    {
				fail <- fail+1
				stop("Error ", tr, " occured. It is likely you are trying to insert duplicate values in the NSV")

			    }
			    total <- total+1

			    if(total%%num_records == 0) cat("CreateNSVTable: Read", total, "entries (ok:", ok, 
				    "/ fail:", fail,")\n")

			} #for(i in 1:nrow(d))
			dbCommit(db$db)
			start<-start+num_records

		    } #while(TRUE)

		    #end loop


		    cat("CreateNSVTable: Read", ok+fail, "entries. Added", ok , "entries.\n")

		}


		dropNSVTable <-  function(db, table) {
		    dbSendQuery(db$db,
			    statement = paste("DROP TABLE ", table, sep='')
			    )
		    dbSendQuery(db$db,statement= paste("DELETE FROM metaData where name='", table,"'",sep=''))
		    dbCommit(db$db)
		    invisible(NULL)
		}

getNSVs <- function(db,  rank=NULL, name=NULL, 
	table="NSV", limit=NULL, random=FALSE, start=1, length=NULL,
	partialMatch=TRUE, removeUnknownSpecies=FALSE, annotation="id") {

    # limit = number of sequences to limit	
    # random = whether the sequences should be random
    # start = start of the chunk eg: 1
    # length = length of the chunk eg: 100 (should be called width?)

    if (!is.null(rank) && rank=="id") partialMatch <- FALSE
    if(is.null(limit)) limit <- "" 
    else limit <- paste(" LIMIT ", paste(limit,collapse=","))

	if(random) limit <- paste(" ORDER BY RANDOM() ", limit)

    if (!is.null(rank)) {    
	fullRank<-BioTools:::.pmatchRank(db,rank)
	#Do this so that the column order appears as [order] since order is a SQL keyword
	fullRankSQL<-paste("classification.[",fullRank,"]",sep="")
    }
    else fullRankSQL <-"-1"

	res <- dbGetQuery(db$db, 
	    statement = paste("SELECT ", table, ".data AS data, classification.id AS id, ", 
		    fullRankSQL ," AS fullRank  FROM sequences ", 
		    " INNER JOIN classification ON classification.id = sequences.id INNER JOIN ", table, " ON ", table, ".id=sequences.id ",
		    BioTools:::.getWhere(db, rank, name, partialMatch,
			    removeUnknownSpecies), limit)
	    )

    if (nrow(res) == 0) stop("No rows found in the database")

    #get metadata about the table
    meta<-as.character(subset(metaGenDB(db),name==table)["annotation"])
    x<-unlist(strsplit(meta,";"))	
    window <-as.numeric(sub("window=","",x[3]))
    overlap <- as.numeric(sub("overlap=","",x[4]))
    word <-as.numeric(sub("word=","",x[5]))
    last_window <-sub("last_window=","",x[6])

    ret <- lapply(res$data,decodeSequence)
    names(ret) <- res$id

    rm(db)

    NSVSet(ret, window=window,overlap=overlap, word=word, 
	    last_window=last_window, rank=rank, name=name)
}


### overwrite getSequences so it can decode NSVs
getSequences <- function(db,  rank=NULL, name=NULL,
	table="sequences", limit=NULL, random=FALSE, start=1, 
	length=NULL, partialMatch=TRUE, removeUnknownSpecies=FALSE, 
	annotation="id") {

    ### FIXME: check metadata table
    type <- metaGenDB(db, table)[,"type"]
    if(length(type) !=1) stop("table does not exist!")    

    if (type=="sequence") getX <- BioTools:::getSequences
    else getX <- getNSVs

    getX(db,  rank, name,
	table, limit, random, start, 
	length, partialMatch, removeUnknownSpecies, 
	annotation)
}
					    
## helper
decodeSequence <- function(sequence) {
    if(length(sequence)==1) unserialize(base64decode(sequence, what="raw"))
    else lapply(sequence, FUN=function(x) unserialize(base64decode(x, what="raw")))
    }




