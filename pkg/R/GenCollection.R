
print.GenCollection <- function(object) {
    cat("Object of class GenCollection")
    ## FIXME: Implement
    #    cat(" of type", object$type)
    cat(" with", nSequences(object), "sequences.\n")
    cat("Collection name:", object$collection, "\n")    
}



createGenDB <- function(db, collection, 
	classification=GenClass16S_Greengenes(),
	type="sequence") {

    ## we use a c_ prefix so SQL is happy
	#first table stores the Classification with org_name as PK
    cl <- paste(paste("'",names(classification),"'", sep=''), "TEXT", collapse=', ')
    cl <- paste(cl, "PRIMARY KEY") ## the lowest rank is the primary key
    
    ## FIXME: NSVs?
    #dat <- "sequence BLOB"
	#first table is called classification and stores the class hierarchy

    try(
		dbSendQuery(db, 
		    statement = paste("CREATE TABLE classification ( " ,
			    cl, ")", sep=''))
	    )
	#second table stores the sequences as BLOB with org_name as PK
	seq <- "org_name TEXT PRIMARY KEY REFERENCES classification(org_name), sequence BLOB "
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
			statement = "INSERT INTO metaData  VALUES ('sequences','sequence','2010 core set')")
		)
	
    openGenCollection(db, collection)
}


openGenCollection <- function(db, collection) {
    col <- list(db=db, collection=collection)   
    class(col) <- "GenCollection"
    col
}

removeGenCollection <- function(db, name) {
	    dbSendQuery(db, 
		    statement = paste("DROP TABLE ", name, sep='')
	    )
}

getClassification <- function(col) {
    cl <- dbListFields(col$db, "classification")
    cl <- head(cl, length(cl)-1L)   ### remove data
    cl
}

nSequences <- function(collection, rank=NULL, name=NULL) {
    dbGetQuery(collection$db, 
	    statement = paste("SELECT COUNT(*) FROM classification", 
	    .getWhere(collection, rank, name)))[1,1]
   
}

getRank <- function(collection, rank=NULL, whereRank=NULL, whereName=NULL) {
    fields <- getClassification(collection)
    cols <- paste("[", fields[.pmatchRank(collection, rank, 
		    numeric=TRUE)],"]", sep='')
    dbGetQuery(collection$db, 
	    statement = paste("SELECT DISTINCT ",cols,
		    " FROM classification ", 
	    .getWhere(collection, whereRank, whereName)))
}
#FIXME: Need to rewrite query to get sequences
getSequences <- function(collection,  rank=NULL, name=NULL, table="sequences",n=-1) {
	table<-sub(" ","",table)
	dbGetQuery(collection$db, 
	    statement = paste("SELECT * FROM ", table ," INNER JOIN classification ON classification.org_name = ", table, ".org_name ", 
		    .getWhere(collection, rank, name))
	    )
	
}


.pmatchRank <- function(col, rank, numeric=FALSE) {
    fields <- dbListFields(col$db, "classification")
    if(is.null(rank)) m <- 1
    else {
	m <- pmatch(tolower(rank), tolower(fields))
	if(any(is.na(m))) stop("Rank not found in collection!")
    }

    if(numeric) m
    else fields[m]
    }

.getWhere <- function(col, rank, name) {
    if(is.null(rank) && is.null(name)) where <- ""
    else where <- paste("WHERE classification.", .pmatchRank(col, rank), 
		" LIKE '", name,"%'", sep='')
    where
}



toNSV <- function(collection, tableName, window=100,
	overlap=0, word=3, last_window=FALSE) {

    ### this might need to much memory. Use SQL LIMIT
	NSV <- "org_name TEXT PRIMARY KEY REFERENCES classification(org_name), NSV BLOB"
	try(
	    dbSendQuery(db, 
		    statement = paste("CREATE TABLE ", tableName , "(",  NSV,  ")", sep='')
	    )
	)
    d <- dbGetQuery(collection$db, 
	    statement = paste("SELECT * FROM sequences"))    
    for(i in 1:nrow(d)) {
	org_name <- d$org_name
	nsv <- counter(d$sequence[i], window, overlap, word, last_window)
	d$sequence[i] <- base64encode(serialize(nsv, NULL))
	#enc <- base64encode(serialize(list(1:10, "A"), NULL, ascii=FALSE))
	#unserialize(base64decode(enc, what="raw"))

	dat <- paste("'",d[i,],"'", sep='', collapse=', ')
	try(
			
		dbSendQuery(collection$db,          
			statement = paste("INSERT INTO ", tableName ," VALUES (", dat, ")", sep=''))
		)
	    
	}
	#insert into meta
	metaAnnotation<- paste("window=",window,"overlap=",overlap,"word=",word,"last_window=",last_window,sep=" ")
	meta<-paste("'",tableName, "','NSV','",metaAnnotation,"'",sep="")
	try(
			
		dbSendQuery(collection$db,          
			statement = paste("INSERT INTO metaData VALUES (", meta,")", sep=''))
		)
	cat("Inserted new table ",tableName," \n")
    
}

decodeSequence <- function(sequence) {
    if(length(sequence)==1) unserialize(base64decode(sequence, what="raw"))
    else lapply(sequence, FUN=function(x) unserialize(base64decode(x, what="raw")))
}

