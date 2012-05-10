## classification hierarchy for 16S
GenClass16S <- function(domain=NA, phylum=NA, class=NA, order=NA, 
	family=NA, genus=NA, species=NA, strain=NA) {

    c(Domain=domain, Phylum=phylum, Class=class,
	    Order=order, Family=family, Genus=genus, Species=species,
	    Strain=strain)
}

### get classification info
getClassification <- function(db) {
    cl <- dbListFields(db$db, "classification")
    cl <- head(cl, length(cl))   ### remove data
    cl
}

getRank <- function(db, rank=NULL, whereRank=NULL, whereName=NULL, all=FALSE) {
    fields <- getClassification(db)
    cols <- paste("[", fields[.pmatchRank(db, rank, 
		    numeric=TRUE)],"]", sep='')

    if(all) distinct <- "" else  distinct <- "DISTINCT"

    dbGetQuery(db$db, 
	    statement = paste("SELECT", distinct, "classification.",cols,
		    " FROM classification ", 
		    .getWhere(db, whereRank, whereName), " ORDER BY ",cols))
}

getHierarchy <- function(db, rank, name)
{
    hierarchy<- GenClass16S_Greengenes()
    rankNum<- which(tolower(names(hierarchy))==tolower(.pmatchRank(db,rank)))
    for(i in 1:(rankNum))
    {
	hierarchy[i]<-getRank(db,rank=names(hierarchy)[i],whereRank=rank,whereName=name)
    }
    #hierarchy[rankNum] <- .pmatchRank(db,rank)
    return(hierarchy)
}

nSequences <- function(db, rank=NULL, name=NULL) {
    dbGetQuery(db$db, 
	    statement = paste("SELECT COUNT(*) FROM classification", 
		    .getWhere(db, rank, name)))[1,1]

}


getSequences <- function(db,  rank=NULL, name=NULL, 
	table="sequences", limit=-1, random=FALSE, start=1, length=NULL) {

    if(limit[1]<0) limit <- "" 
    else limit <- paste(" LIMIT ",paste(limit,collapse=","))

	if(random)  
	limit <- paste(" ORDER BY RANDOM() ",limit)
    #make length SQL compatible
    if (is.null(length)) lengthFilter= "data"
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
			.getWhere(db, rank, name), limit)
		)
    }
    else {	    
	res <- dbGetQuery(db$db, 
		statement = paste("SELECT ",lengthFilter ,"  AS data, classification.id AS id, ", fullRankSQL ," AS fullRank  FROM ", table ,
			" INNER JOIN classification ON classification.id = ",
			table, ".id ", 
			.getWhere(db, rank, name), limit)
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


