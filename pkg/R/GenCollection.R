
print.GenCollection <- function(object) {
    cat("Object of class GenCollection")
    ## FIXME: Implement
    #    cat(" of type", object$type)
    cat(" with", nSequences(object), "sequences.\n")
    cat("Collection name:", object$collection, "\n")    
}


createGenCollection <- function(db, collection, 
	classification=GenClass16S_Greengenes(),
	type="sequence") {

    ## we use a c_ prefix so SQL is happy
    cl <- paste(paste("'",names(classification),"'", sep=''), "TEXT", collapse=', ')
    cl <- paste(cl, "PRIMARY KEY") ## the lowest rank is the primary key
    
    ## FIXME: NSVs?
    dat <- "sequence BLOB"

    try(
	    dbSendQuery(db, 
		    statement = paste("CREATE TABLE ", collection, '(',
			    cl, ', ', dat, ')', sep=''))
	    )

    openGenCollection(db, collection)
}

openGenCollection <- function(db, collection) {
    col <- list(db=db, collection=collection)   
    class(col) <- "GenCollection"
    col
}





getClassification <- function(col) {
    cl <- dbListFields(col$db, col$collection)
    cl <- head(cl, length(cl)-1L)   ### remove data
    cl
}



nSequences <- function(collection, rank=NULL, name=NULL) {
    dbGetQuery(collection$db, 
	    statement = paste("SELECT COUNT(*) FROM", collection$collection, 
	    .getWhere(collection, rank, name)))[1,1]
   
}

getRank <- function(collection, rank=NULL, whereRank=NULL, whereName=NULL) {
    fields <- getClassification(collection)
    cols <- paste("[", fields[.pmatchRank(collection, rank, 
		    numeric=TRUE)],"]", sep='')
    dbGetQuery(collection$db, 
	    statement = paste("SELECT DISTINCT ",cols,
		    " FROM", collection$collection, 
	    .getWhere(collection, whereRank, whereName)))
}


.pmatchRank <- function(col, rank, numeric=FALSE) {
    fields <- dbListFields(col$db, col$collection)
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
    else where <- paste("WHERE ", .pmatchRank(col, rank), 
		" LIKE '", name,"%'", sep='')
    where
}

getSequences <- function(collection, rank=NULL, name=NULL, n=-1) {
    dbGetQuery(collection$db, 
	    statement = paste("SELECT * FROM ", collection$collection, " ", 
		    .getWhere(collection, rank, name))
	    )
}

toNSV <- function(from, to, window=100,
	overlap=0, word=3, last_window=FALSE) {

    ### this might need to much memory. Use SQL LIMIT
    d <- dbGetQuery(from$db, 
	    statement = paste("SELECT * FROM ", from$collection))
    
    for(i in 1:nrow(d)) {
	nsv <- counter(d$sequence[i], window, overlap, word, last_window)
	d$sequence[i] <- base64encode(serialize(nsv, NULL))
	#enc <- base64encode(serialize(list(1:10, "A"), NULL, ascii=FALSE))
	#unserialize(base64decode(enc, what="raw"))

	dat <- paste("'",d[i,],"'", sep='', collapse=', ')
	dbSendQuery(to$db,          
		statement = paste("INSERT INTO ",
			to$collection, " VALUES(", 
			dat, ")", sep=''))

    }    
}
    
decodeSequence <- function(sequence) {
    if(length(sequence)==1) unserialize(base64decode(sequence, what="raw"))
    else lapply(sequence, FUN=function(x) unserialize(base64decode(x, what="raw")))
}

