### print an NSVSet object
print.NSVSet <- function(x, ...) {
    cat("Object of class NSVSet for", length(x), "sequences")
    cat(" (", log(length(colnames(x[[1]])),4) ,"-mers)\n", sep="")
    cat("Number of segments (table with counts):") 
    print(table(sapply(x, nrow)))
}

### subset
`[.NSVSet` <- function(x, i, j, ..., drop = TRUE) { 
    r <- unclass(x)
    if(!missing(i)) r <- r[i] 
    if(!missing(j)) r <- lapply(r, "[", j, ,drop=FALSE) 
    class(r) <- "NSVSet"
    r
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
    class(s) <- "NSVSet"
    s
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
    num_records<-100 # number of records at a time
	if (is.null(limit))
		limit <- 500000
    while(TRUE && total <limit)
    {
	if (is.null(rank) && is.null(name))
	    d <- dbGetQuery(db$db, statement = 
		paste("SELECT * FROM sequences ", if (removeUnknownSpecies) " INNER JOIN classification c ON c.id=sequences.id WHERE c.[Species] NOT LIKE 'Unknown%'" ,"LIMIT ", start,",
			",num_records,sep=""))    
	else {
		
		#statement<-	paste("SELECT sequences.id, sequences.data FROM sequences INNER JOIN classification ON sequences.id=classification.id WHERE classification.",
		#	.pmatchRank(db,rank)," LIKE '", 
		#	name,"%'", if(removeUnknownSpecies) " AND classification.[Species] NOT LIKE 'Unknown%' ", " LIMIT ", start,", ",num_records, sep="")
		statement<-	paste("SELECT sequences.id, sequences.data FROM sequences INNER JOIN classification ON sequences.id=classification.id ", .getWhere(db, rank, name) ,
			if(removeUnknownSpecies) " AND classification.[Species] NOT LIKE 'Unknown%' ", " LIMIT ", start,", ",num_records, sep="")
	
		
    d <- dbGetQuery(db$db, statement = 
		paste("SELECT sequences.id, sequences.data FROM sequences INNER JOIN classification ON sequences.id=classification.id ", .getWhere(db, rank, name) ,
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
    invisible()
}

### helper
decodeSequence <- function(sequence) {
    if(length(sequence)==1) unserialize(base64decode(sequence, what="raw"))
    else lapply(sequence, FUN=function(x) unserialize(base64decode(x, what="raw")))
    }




