## read fasta files and add them to a DB

addSequences <- function(db, file, metaDataReader=GreengenesMetaDataReader, 
	verbose=FALSE) {

    ok <- 0
    fail <- 0
    total <- 0

    dbBeginTransaction(db$db)
    #start
	f <- read.DNAStringSet(file)
	for(i in 1:length(f)) {
		annot<- names(f)[i]
		cl <- metaDataReader(annot)
		org_name<-cl[length(cl)]
		cl <- paste("'",cl,"'", sep='', collapse=', ') 
		#print(cl)
		dat<- f[[i]]
		dat<- tolower(as.character(dat[1:length(dat)]))
		tr <- try(dbSendQuery(db$db,          
			statement = paste("INSERT INTO classification VALUES(", 
				cl,  ")", sep='')), silent=FALSE)
		
		tr <- try(dbSendQuery(db$db,          
			statement = paste("INSERT INTO sequences VALUES('", 
				org_name, "','", dat,  "')", sep='')), 
				silent=FALSE)
		if(!is(tr, "try-error")) ok <- ok+1
		else fail <- fail+1
		total <- total+1
		if(verbose)
		{	
			if(total%%100 == 0) cat("Read", total, "entries (ok:", ok, 
				"/ fail:", fail,")\n")
		}
	}


    dbCommit(db$db)

    #close(f)
    
    cat("Read", ok+fail, "entries. Added", ok , "entries.\n")

}



