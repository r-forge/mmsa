### Sequences are DNAStringSet

### plot as a Sequence Logo
plot.DNAMultipleAlignment <- function(x, start=1, end=NULL, ic.scale=FALSE, ...) plot(as(x, "DNAStringSet"), start, end, ic.scale, ...) 


plot.DNAStringSet <- function(x, start=1, end=NULL, ic.scale=FALSE, ...) 
{
    cm <- consensusMatrix(x,baseOnly=TRUE)
    if (is.null(end))
	end <- ncol(cm)
    cm <- cm[1:4,start:end]
    for(i in 1:ncol(cm))
	cm[,i] <- cm[,i]/colSums(cm)[i]
    pwm <- makePWM(cm)
    seqLogo::seqLogo(pwm, ic.scale,  ...)
}

### print, [, etc. is handeled by Biostrings

### GenDB
nSequences <- function(db, rank=NULL, name=NULL) {
    dbGetQuery(db$db, 
	    statement = paste("SELECT COUNT(*) FROM classification", 
		    .getWhere(db, rank, name)))[1,1]

}


getSequences <- function(db,  rank=NULL, name=NULL, 
	table="sequences", limit=NULL, random=FALSE, start=1, length=NULL,
	partialMatch=TRUE) {

	# limit = number of sequences to limit	
	# random = whether the sequences should be random
	# start = start of the chunk eg: 1
	# length = length of the chunk eg: 100

    if(is.null(limit)) limit <- "" 
    else limit <- paste(" LIMIT ",paste(limit,collapse=","))

	if(random)  
		limit <- paste(" ORDER BY RANDOM() ",limit)
    #get chunks of sequences, important for clustering
	#make length SQL compatible
    if (is.null(length)) 
		lengthFilter= "data"
    else
		lengthFilter = paste("substr(data,",start,",",length,")",sep="")

    if (!is.null(rank)) {    
		fullRank<-.pmatchRank(db,rank)
		#Do this so that the column order appears as [order] since order is a SQL keyword
		fullRankSQL<-paste("classification.[",fullRank,"]",sep="")
    }
    else
	fullRankSQL <-"-1"

    #different route for NSV segments
    #get Sequences in memory and convert to NSV using in-memory function createNSVSet
    if(table!="sequences" && !is.null(length)) {
	res <- dbGetQuery(db$db, 
		statement = paste("SELECT ", lengthFilter ," AS data, classification.id AS id, ", fullRankSQL ," AS fullRank  FROM sequences ", 
			" INNER JOIN classification ON classification.id = sequences.id ",
			.getWhere(db, rank, name, partialMatch), limit)
		)
    }
    else {	    
	res <- dbGetQuery(db$db, 
		statement = paste("SELECT ",lengthFilter ,"  AS data, classification.id AS id, ", fullRankSQL ," AS fullRank  FROM ", table ,
			" INNER JOIN classification ON classification.id = ",
			table, ".id ", 
			.getWhere(db, rank, name, partialMatch), limit)
		)
    }

    if (nrow(res) == 0) stop("No rows found in the database")
    if (table !="sequences") {
	#get metadata about the table
		meta<-as.character(subset(metaGenDB(db),name==table)["annotation"])
		x<-unlist(strsplit(meta,";"))	
		window <-as.numeric(sub("window=","",x[3]))
		overlap <- as.numeric(sub("overlap=","",x[4]))
		word <-as.numeric(sub("word=","",x[5]))
		last_window <-sub("last_window=","",x[6])
		#need to convert DNA sequences into NSV in-memory
		if(!is.null(length))
	    	ret <- createNSVSet(res$data, window=window, overlap=overlap, word=word, last_window=last_window)   
		else    
	    	ret <- lapply(res$data,decodeSequence)
		attr(ret,"window")<- window
		attr(ret,"overlap") <- overlap
		attr(ret,"word")<- word
		attr(ret,"last_window")<- last_window
		if(!is.null(rank)) {
	    	attr(ret,"rank")<-fullRank
	    	#this returns the values of the fullRank i.e. rankName from the db
	    	attr(ret,"name")<-res$fullRank
		}
		attr(ret,"id") <- res$id
		names(ret) <- res$id
		class(ret) <- "NSVSet"
    } else {
	ret <- DNAStringSet(res$data)
	names(ret) <- res$id
	if(!is.null(rank)){
	    attr(ret,"rank")<-fullRank
	    attr(ret,"name")<-res$fullRank
	}
    }
    ret
}


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


