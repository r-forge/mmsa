
createGenDB <- function(dbName, classification=GenClass16S_Greengenes(),
	 drv=NULL, ...) {

    if(file.exists(dbName)) 
	stop("GenDB already exist. Use openGenDB.\n")
    
    
    if(is.null(drv)) drv<-dbDriver("SQLite");
    db<-dbConnect(drv, dbname = dbName, ...);


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

    db <- list(db=db, dbName=dbName, drv=drv)
    class(db) <- "GenDB"
    db
}

openGenDB <- function(dbName, drv=NULL, ...) {
    if(!file.exists(dbName)) 
	stop("GenDB does not exist. Use createGenDB first.\n")
    
    if(is.null(drv)) drv<-dbDriver("SQLite");

    db<-dbConnect(drv, dbname = dbName, ...);

    db <- list(db=db, dbName=dbName, drv=drv)
    class(db) <- "GenDB"
    db
}

closeGenDB <- function(db) {
    dbCommit(db$db)
    ret <- dbDisconnect(db$db)
    return(invisible(ret))
}

reopenGenDB <- function(db, ...) {
    db_new<-dbConnect(db$drv, dbname = db$dbName, ...);
    db <- list(
	    db=dbConnect(db$drv, dbname = db$dbName), 
	    dbName=db$dbName, drv=db$drv)
    class(db) <- "GenDB"
    db
}

#reopenGenDB <- function(db, ...) {
#    db_new<-dbConnect(db$drv, dbname = db$dbName, ...);
#    db <- list(
#	    db=db_new, 
#	    dbName=db$dbName, drv=db$drv)
#    class(db_new) <- "GenDB"
#    db_new
#}

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
    cat(paste(dbListTables(x$db), collapse=", "), "\n")
}

randomizeGenDB <- function(db)
{
    tables <- listGenDB(db)
    for(i in 1:length(tables))
    {
        dbSendQuery(db$db,statement= paste("CREATE TABLE temp AS SELECT * FROM ",tables[i]," ORDER BY random()",sep=''))
        dbSendQuery(db$db,statement= paste("DROP TABLE ",tables[i],sep=''))
        dbSendQuery(db$db,statement= paste("ALTER TABLE temp RENAME TO  ",tables[i],sep=''))
    }
    dbCommit(db$db)
}

