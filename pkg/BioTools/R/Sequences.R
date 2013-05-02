#######################################################################
# BioTools - Interfaces to several sequence alignment and 
# classification tools
# Copyright (C) 2012 Michael Hahsler and Anurag Nagar
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.


### Sequences are DNAStringSet

### GenDB
nSequences <- function(db, rank=NULL, name=NULL, table="sequences") {
    dbGetQuery(db$db, 
	    statement = paste("SELECT COUNT(*) FROM ",table," t INNER JOIN classification ON t.id=classification.id ", 
		    .getWhere(db, rank, name)))[1,1]

}


getSequences <- function(db,  rank=NULL, name=NULL, 
	table="sequences", limit=NULL, random=FALSE, start=1, length=NULL,
	partialMatch=TRUE, removeUnknownSpecies=FALSE, annotation="id") {

    # limit = number of sequences to limit	
    # random = whether the sequences should be random
    # start = start of the chunk eg: 1
    # length = length of the chunk eg: 100 (should be called width?)

    if (!is.null(rank))
	if(rank=="id") partialMatch <- FALSE
    if(is.null(limit)) limit <- "" 
    else limit <- paste(" LIMIT ",paste(limit,collapse=","))

	if(random)  
	limit <- paste(" ORDER BY RANDOM() ",limit)
    #get chunks of sequences, important for clustering
    #make length SQL compatible
    if (is.null(length)) 
	lengthFilter= "data"
    else
	lengthFilter = paste("substr(sequences.data,",start,",",length,")",
	    sep="")

    if (!is.null(rank)) {    
	fullRank<-.pmatchRank(db,rank)
	#Do this so that the column order appears as [order] since order is a SQL keyword
	fullRankSQL<-paste("classification.[",fullRank,"]",sep="")
    }
    else fullRankSQL <-"-1"

	res <- dbGetQuery(db$db, 
	    statement = paste("SELECT ",lengthFilter ,"  AS data, classification.id AS id, ", fullRankSQL ," AS fullRank  FROM ", table ,
		    " INNER JOIN classification ON classification.id = ",
		    table, ".id ", 
		    .getWhere(db, rank, name, partialMatch,removeUnknownSpecies), limit)
	    )
	if (!is.null(rank) && rank=="id")
	{
		 if (!is.null(name))
			res<-res[match(name,res$id),]	
	}
	if (nrow(res) == 0) stop("No rows found in the database")
    ret <- DNAStringSet(res$data)
    if (annotation=="id")
	names(ret) <- res$id
    else if (annotation=="rdp")
    {
	h <- getHierarchy(db, rank="id", name=res$id, partialMatch=FALSE)[,1:6]
	if (class(h)=="matrix")
	    hierarchy <- apply(h,MARGIN=1,FUN=function(x) paste(x,collapse=";"))
	hierarchy <- gsub(";unknown","",hierarchy)
	hierarchy <- gsub(" \\(class\\)","",hierarchy)
	hierarchy <- paste("Root",hierarchy,sep=";")
	hierarchy <- paste(res$id,hierarchy)	
	names(ret) <- hierarchy
    }
    if(!is.null(rank)){
	attr(ret,"rank")<-fullRank
	attr(ret,"name")<-res$fullRank
    }

    ret
}


## read fasta files and add them to a DB

addSequences <- function(db, file, metaDataReader=GreengenesMetaDataReader, 
	verbose=FALSE) {

    if(!file.exists(file)) stop("File does not exist!")
    if(file.info(file)$isdir) {
	cat("Found directory. Adding whole directory.\n")
	file <- list.files(file, full.names=TRUE, 
		recursive=TRUE)
    }


    ok <- 0
    fail <- 0
    total <- 0

    dbBeginTransaction(db$db)
    #start
	f <- readDNAStringSet(file)
	for(i in 1:length(f)) {
		annot<- names(f)[i]
		cl <- metaDataReader(annot)
		org_name<-cl[length(cl)]
		cl <- paste("'",cl,"'", sep='', collapse=', ') 
		dat<- f[[i]]
		dat<- tolower(as.character(dat[1:length(dat)]))
		try(dbSendQuery(db$db,          
			statement = paste("INSERT INTO classification VALUES(", 
				cl,  ")", sep='')), silent=TRUE)
		
		tr <- try(dbSendQuery(db$db,          
			statement = paste("INSERT INTO sequences VALUES('", 
				org_name, "','", dat,  "')", sep='')), 
				silent=TRUE)
		
		if(is(tr, "try-error")) { 
		    if(verbose) cat("Adding", org_name, "failed -",
			    attr(tr, "condition")$message, "\n")
		    fail <- fail+1
		}else ok <- ok+1
		
		
		total <- total+1
		if(verbose)
		{	
			if(total%%100 == 0) cat("Read", total, "sequences (ok:", ok, 
				"/ fail:", fail,")\n")
		}
	}

    dbCommit(db$db)
    #close(f)
    cat("Read", ok+fail, "sequences. Added", ok , "sequences.\n")
}


