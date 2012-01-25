### convert to NSVs

createNSVTable <- function(db, tableName, whereRank=NULL,whereName=NULL, window=100,
	overlap=0, word=3, last_window=FALSE) {

    if(length(grep(" ", tableName))) stop("tableName cannot contain spaces!")

    ### this might need to much memory. Use SQL LIMIT
    NSV <- "id TEXT PRIMARY KEY REFERENCES classification(id), data BLOB"
    try(
	    dbSendQuery(db$db, 
		    statement = paste("CREATE TABLE ", tableName ,
			    "(",  NSV,  ")", sep='')
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
		if (is.null(whereRank) && is.null(whereName))
			d <- dbGetQuery(db$db, statement = paste("SELECT * FROM sequences LIMIT ",start,",",num_records,sep=""))    
		else
			d <- dbGetQuery(db$db, statement = paste("SELECT sequences.id, sequences.data FROM sequences INNER JOIN classification ON sequences.id=classification.id WHERE classification.",
					.pmatchRank(db,whereRank)," LIKE '", whereName,"%' LIMIT ",start,",",num_records,sep="" ))    
		
		if (nrow(d)==0) break;

		dbBeginTransaction(db$db)
		for(i in 1:nrow(d))  {
	
		#make NSV
	
		nsv <- .counter(d$data[i], window, overlap, word, last_window)
		d$data[i] <- base64encode(serialize(nsv, NULL))
		#end make NSV
		#org_name<-d$org_name_
	
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
				statement = paste("INSERT INTO ", tableName," VALUES (", dat, ")", sep='')
				#print(statement)
			}
		total <- total+1

		if(total%%num_records == 0) cat("Read", total, "entries (ok:", ok, 
			"/ fail:", fail,")\n")

		} #for(i in 1:nrow(d))
		dbCommit(db$db)
		start<-start+num_records
		
	} #while(TRUE)
		
	#end loop
    
    
    cat("Read", ok+fail, "entries. Added", ok , "entries.\n")

}

#end


dropNSVTable <-  function(db, tableName) {
    dbSendQuery(db$db,
	    statement = paste("DROP TABLE ", tableName, sep='')
	    )
    invisible()
}

### helper
decodeSequence <- function(sequence) {
    if(length(sequence)==1) unserialize(base64decode(sequence, what="raw"))
    else lapply(sequence, FUN=function(x) unserialize(base64decode(x, what="raw")))
    }




