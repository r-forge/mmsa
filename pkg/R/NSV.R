### convert to NSVs

createNSVTable <- function(db, tableName, window=100,
	overlap=0, word=3, last_window=FALSE) {

    if(length(grep(" ", tableName))) stop("tableName cannot contain spaces!")

    ### this might need to much memory. Use SQL LIMIT
    NSV <- "org_name TEXT PRIMARY KEY REFERENCES classification(org_name), NSV BLOB"
    try(
	    dbSendQuery(db$db, 
		    statement = paste("CREATE TABLE ", tableName ,
			    "(",  NSV,  ")", sep='')
		    )
	    )
    
    
    d <- dbGetQuery(db$db, 
	    statement = paste("SELECT * FROM sequences"))    
    
    
    for(i in 1:nrow(d)) {
	org_name <- d$org_name
	nsv <- counter(d$sequence[i], window, overlap, word, last_window)
	d$sequence[i] <- base64encode(serialize(nsv, NULL))
	
	#enc <- base64encode(serialize(list(1:10, "A"), NULL, ascii=FALSE))
	#unserialize(base64decode(enc, what="raw"))

	dat <- paste("'",d[i,],"'", sep='', collapse=', ')
	
	try(
		dbSendQuery(db$db,          
			statement = paste("INSERT INTO ", tableName,
				" VALUES (", dat, ")", sep=''))
		)

    }
    
    #insert into meta
    metaAnnotation<- paste("window=", window, ", overlap=", overlap,
	    ", word=", word, ", last_window=", last_window, sep="")
    meta<-paste("'", tableName, "','NSV','", metaAnnotation, "'", sep="")
    
    try(
	    dbSendQuery(db$db,          
		    statement = paste("INSERT INTO metaData VALUES (", meta,")", sep=''))
	    )
    
    #cat("Inserted new table ", tableName," \n")
    invisible()

}

createNSVTable_large <- function(db, tableName, window=100,
	overlap=0, word=3, last_window=FALSE) {

    if(length(grep(" ", tableName))) stop("tableName cannot contain spaces!")

    ### this might need to much memory. Use SQL LIMIT
    NSV <- "org_name TEXT PRIMARY KEY REFERENCES classification(org_name), NSV BLOB"
    try(
	    dbSendQuery(db$db, 
		    statement = paste("CREATE TABLE ", tableName ,
			    "(",  NSV,  ")", sep='')
		    )
	    )
    ok <- 0
    fail <- 0
    total <- 0

	d <- dbGetQuery(db$db, statement = paste("SELECT * FROM sequences"))    

    dbBeginTransaction(db$db)
   
	
	for(i in 1:nrow(d))  {
	
	#make NSV
	
	nsv <- counter(d$sequence[i], window, overlap, word, last_window)
	d$sequence[i] <- base64encode(serialize(nsv, NULL))
	#end make NSV
	org_name<-d$org_name
	
	## Insert into DB
	
	dat <- paste("'",d[i,],"'", sep='', collapse=', ')
	
	tr<- try(
		dbSendQuery(db$db,          
			statement = paste("INSERT INTO ", tableName,
				" VALUES (", dat, ")", sep=''))
		)

    if(!is(tr, "try-error")) ok <- ok+1
	else fail <- fail+1
	total <- total+1

	if(total%%100 == 0) cat("Read", total, "entries (ok:", ok, 
		"/ fail:", fail,")\n")

    } #for(i in 1:nrow(d))
    dbCommit(db$db)

    
    
    cat("Read", ok+fail, "entries. Added", ok , "entries.\n")

}

#end


dropNSVTable <-  function(db, name) {
    dbSendQuery(db$db,
	    statement = paste("DROP TABLE ", name, sep='')
	    )
    invisible()
}

### helper
decodeSequence <- function(sequence) {
    if(length(sequence)==1) unserialize(base64decode(sequence, what="raw"))
    else lapply(sequence, FUN=function(x) unserialize(base64decode(x, what="raw")))
    }




