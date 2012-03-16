
createGenDB <- function(dbName, classification=GenClass16S_Greengenes(),
	 drv=NULL) {

    if(file.exists(dbName)) 
	stop("GenDB already exist. Use openGenDB.\n")
    
    
    if(is.null(drv)) drv<-dbDriver("SQLite");
    db<-dbConnect(drv, dbname = dbName);


    #first table stores the Classification with org_name as PK
    cl <- paste(paste("'",names(classification),"'", sep=''), "TEXT", 
	    collapse=', ')
    cl <- paste(cl, "PRIMARY KEY") ## the lowest rank is the primary key

    #first table is called classification and stores the class hierarchy
    try(
	    dbSendQuery(db, 
		    statement = paste("CREATE TABLE classification ( " ,
			    cl, ")", sep=''))
	    )

    #second table stores the sequences as BLOB with org_name as PK
    seq <- "id TEXT PRIMARY KEY REFERENCES classification(id), data BLOB "
    try(
	    dbSendQuery(db, 
		    statement = paste("CREATE TABLE sequences (",
			    seq,  ")", sep=''))
	    )


    #third table stores the meta data
    meta<-"name TEXT, type TEXT, annotation TEXT"
    try(
	    dbSendQuery(db, 
		    statement = paste("CREATE TABLE metaData (",
			    meta,  ")", sep=''))
	    )

    #insert data into meta
    try(

	    dbSendQuery(db,          
		    statement = paste("INSERT INTO metaData  VALUES ('sequences', 'sequence', '')", sep="")
		    )
	    )

    db <- list(db=db, dbName=dbName)
    class(db) <- "GenDB"
    db
}

openGenDB <- function(dbName, drv=NULL) {
    if(!file.exists(dbName)) 
	stop("GenDB does not exist. Use createGenDB first.\n")
    
    if(is.null(drv)) drv<-dbDriver("SQLite");

    db<-dbConnect(drv, dbname = dbName);

    db <- list(db=db, dbName=dbName)
    class(db) <- "GenDB"
    db
}


closeGenDB <- function(db) {
    dbCommit(db$db)
    ret <- dbDisconnect(db$db)
    return(invisible(ret))
}

listGenDB <- function(db) {
    dbListTables(db$db)
}

metaGenDB <- function(db) {
	dbGetQuery(db$db,"SELECT * from metaData")
}

print.GenDB <- function(x, ...) {
    cat("Object of class GenDB with", nSequences(x), "sequences\n")
    cat("DB File:", x$dbName, "\n")
    cat("Tables: ")
    cat(paste(dbListTables(x$db), ", ", sep=""), "\n")
}


### get classification info

getClassification <- function(db) {
    cl <- dbListFields(db$db, "classification")
    cl <- head(cl, length(cl)-1L)   ### remove data
    cl
}

getRank <- function(db, rank=NULL, whereRank=NULL, whereName=NULL) {
    fields <- getClassification(db)
    cols <- paste("[", fields[.pmatchRank(db, rank, 
		    numeric=TRUE)],"]", sep='')
    dbGetQuery(db$db, 
	    statement = paste("SELECT DISTINCT classification.",cols,
		    " FROM classification ", 
		    .getWhere(db, whereRank, whereName)))
}

getHierarchy <- function(db, rank=NULL, whereRank=NULL, whereName=NULL) {
    fields <- getClassification(db)
    cols <- paste("[", fields[.pmatchRank(db, rank, 
		    numeric=TRUE)],"]", sep='')
    dbGetQuery(db$db, 
	    statement = paste("SELECT  classification.",cols,
		    " FROM classification ", 
		    .getWhere(db, whereRank, whereName)))
}

nSequences <- function(db, rank=NULL, name=NULL) {
    dbGetQuery(db$db, 
	    statement = paste("SELECT COUNT(*) FROM classification", 
		    .getWhere(db, rank, name)))[1,1]

}


getSequences <- function(db,  rank=NULL, name=NULL, table="sequences", limit=-1, random=FALSE) {

    if(limit[1]<0) limit <- "" 
    else limit <- paste(" LIMIT ",paste(limit,collapse=","))
    
	if(random)  
    	limit <- paste(" ORDER BY RANDOM() ",limit)

	if (!is.null(rank)) {    
		fullRank<-.pmatchRank(db,rank)
		#Do this so that the column order appears as [order] since order is a SQL keyword
		fullRankSQL<-paste("[",fullRank,"]",sep="")
	}
	else
		fullRankSQL <-"-1"
	res <- dbGetQuery(db$db, 
		statement = paste("SELECT data, classification.id AS id, ", fullRankSQL ," AS fullRank  FROM ", table ,
		" INNER JOIN classification ON classification.id = ",
		table, ".id ", 
		.getWhere(db, rank, name), limit)
	    )

    if (table !="sequences") {
	ret <- lapply(res$data,decodeSequence)		
	if(!is.null(rank)) {
	    attr(ret,"rank")<-fullRank
	    attr(ret,"name")<-res$fullRank
	}
	attr(ret,"id")<-res$id
	class(ret) <- "NSVSet"
    }else{

	ret <- DNAStringSet(res$data)
	names(ret) <- res$id
	if(!is.null(rank)){
	    attr(ret,"rank")<-fullRank
	    attr(ret,"name")<-res$fullRank
	}
    }
    ret
}


### helper

.pmatchRank <- function(col, rank, numeric=FALSE) {
    fields <- dbListFields(col$db, "classification")
    if(is.null(rank)) m <- 1
    else {
	m <- pmatch(tolower(rank), tolower(fields))
	if(any(is.na(m))) stop("Rank not found!")
    }

    if(numeric) m
    else fields[m]
    }

.getWhere <- function(col, rank, name) {
    if(is.null(rank) && is.null(name)) where <- ""
    else where <- paste("WHERE classification.'", .pmatchRank(col, rank), 
		"' LIKE '", name,"%'", sep='')
    where
}

