### print an NSVSet object
print.NSVSet <- function(x, ...) {
    cat("Object of class NSVSet for", length(x), "sequences")
    cat(" (", log(length(colnames(x[[1]])),4) ,"-mers)\n", sep="")
    cat("number of segments (range):", range(sapply(x, nrow)), "\n")
}

### convert to NSVs
createNSVSet <- function(x, window=100, overlap=0, word=3, 
	last_window=FALSE) {
    s <- lapply(x, .counter, window, overlap, word, last_window)
    class(s) <- "NSVSet"
    s
}


createNSVTable <- function(db, tableName, 
	rank=NULL, name=NULL, window=100,
	overlap=0, word=3, last_window=FALSE) {

    if(length(grep(" ", tableName))) stop("tableName cannot contain spaces!")
    if (length(which(tableName==listGenDB(db))) > 0)
	stop("A table with this name already exists in the db")

    NSV <- "id TEXT PRIMARY KEY REFERENCES classification(id), data BLOB"
    try(
	    dbSendQuery(db$db, 
		    statement = paste("CREATE TABLE ", tableName ,
			    "(",  NSV,  ")", sep='')
		    )
	    )
    # insert into metadata
    meta<-paste("'", tableName , "','NSV','rank=", rank, ";name=",
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
    while(TRUE)
    {
	if (is.null(rank) && is.null(name))
	    d <- dbGetQuery(db$db, statement = 
		paste("SELECT * FROM sequences LIMIT ", start,",
			",num_records,sep=""))    
	else
	    d <- dbGetQuery(db$db, statement = 
		paste("SELECT sequences.id, sequences.data FROM sequences INNER JOIN classification ON sequences.id=classification.id WHERE classification.",
			.pmatchRank(db,rank)," LIKE '", 
			name,"%' LIMIT ", start,", ",num_records, sep=""))    

	if (nrow(d)==0) break;

	dbBeginTransaction(db$db)
	for(i in 1:nrow(d))  {

	    #make NSV

	    nsv <- .counter(d$data[i], window, overlap, word, last_window)
	    d$data[i] <- base64encode(serialize(nsv, NULL))
	    #end make NSV
	    ## Insert into DB
	    dat <- paste("'",d[i,],"'", sep='', collapse=', ')

	    tr<- try(
		    dbSendQuery(db$db,          
			    statement = paste("INSERT INTO ", tableName,
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

#end


dropNSVTable <-  function(db, tableName) {
    dbSendQuery(db$db,
	    statement = paste("DROP TABLE ", tableName, sep='')
	    )
    dbCommit(db$db)
    invisible()
}

### helper
decodeSequence <- function(sequence) {
    if(length(sequence)==1) unserialize(base64decode(sequence, what="raw"))
    else lapply(sequence, FUN=function(x) unserialize(base64decode(x, what="raw")))
    }




