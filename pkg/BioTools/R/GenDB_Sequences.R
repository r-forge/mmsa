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
nSequences <- function(db, rank=NULL, name=NULL, table="classification") {
    dbGetQuery(db$db, 
	    statement = paste("SELECT COUNT(*) FROM ",table,
            " t INNER JOIN classification ON t.id=classification.id ", 
		    .getWhere(db, rank, name)))[1,1]

}

getSequences <- function(db,  rank=NULL, name=NULL, 
	table="sequences", limit=NULL, random=FALSE, start=NULL, length=NULL,
	removeUnknownSpecies=FALSE, annotation=Annotation_Id) {

    # limit = number of sequences to limit	
    # random = whether the sequences should be random
    # start = start of the chunk eg: 1
    # length = length of the chunk eg: 100 (should be called width?)

  if(!is.null(limit)) limitSQL <- paste(" LIMIT ",paste(limit,collapse=","))
  else limitSQL <- "" 
  
  if(random) limitSQL <- paste(" ORDER BY RANDOM() ", limitSQL)
  
  #get chunks of sequences, important for clustering
  #make length SQL compatible
  if (is.null(start) && is.null(length)) lengthFilter= "data"
  else {
    if(is.null(length)) length <- .Machine$integer.max
    if(is.null(start)) start <- 1L
    lengthFilter = paste("SUBSTR(sequences.data,",
                         as.integer(start), ",",
                         as.integer(length), ")", sep="")
  }

  if (!is.null(rank)) {    
    fullRank <- .pmatchRank(db,rank)
    #Do this so that the column order appears as 'order' since ORDER is a SQL keyword
    fullRankSQL<-paste("classification.'",fullRank,"'",sep="")
	}else fullRankSQL <-"-1"

  statement <- paste("SELECT ", lengthFilter, " AS data, classification.id AS id, ", 
                     fullRankSQL ," AS fullRank  FROM ", table ,
                     " INNER JOIN classification ON classification.id = ",
                     table, ".id ", 
                     .getWhere(db, rank, name, removeUnknownSpecies), 
                     limitSQL, sep='')
  
	res <- dbGetQuery(db$db, statement = statement)

  if (!is.null(rank) && rank=="id" && nrow(res)==length(name)) {
		if (!is.null(name)) res<-res[match(name,res$id),]	
	}
	
  if (nrow(res) == 0) stop("No rows found in the database")
  ret <- DNAStringSet(res$data)
  
  if(identical(annotation, Annotation_Id)) names(ret) <- res$id
  else { ### this is slow
    h <- getHierarchy(db, rank="id", name=res$id, 
                      drop=FALSE)
    names(ret) <- annotation(h, decode=FALSE)
  }
    
  if(!is.null(rank)){
    attr(ret,"rank")<-fullRank
    attr(ret,"name")<-res$fullRank
  }
  
  ret
}

## read fasta files or a DNAStringSet and add them to a DB
addSequences <- function(db, sequences, table="sequences",
                         annotation=Annotation_Greengenes, verbose=FALSE) {
  
  src <- deparse(substitute(sequences))
  
  if(!is(sequences, "DNAStringSet")) {
    if(!file.exists(sequences)) stop("File does not exist!")
    
    if(file.info(sequences)$isdir) {
      cat("Found directory. Adding whole directory.\n")
      sequences <- list.files(sequences, full.names=TRUE, 
                         recursive=TRUE)
    }
    
    sequences <- readDNAStringSet(sequences)
  }
  
  ### does table exist?
  if(!(table %in% listGenDB(db))) 
    .createdatatableGenDB(db, table, "sequence", annotation=paste("source:", src))
  
  ok <- 0
  fail <- 0
  total <- 0
  
  dbBeginTransaction(db$db)
  #start
  for(i in 1:length(sequences)) {
    cl <- annotation(names(sequences)[i], decode=TRUE)
    id <- cl["Id"]
    cl <- paste("'",cl,"'", sep='', collapse=', ') 
    dat<- as.character(sequences[[i]])
    
    try(dbSendQuery(db$db,          
                    statement = paste("INSERT INTO classification VALUES(", 
                                      cl,  ")", sep='')), silent=TRUE)
    
    tr <- try(dbSendQuery(db$db,          
                          statement = paste("INSERT INTO ", table, " VALUES('", 
                                            id, "','", dat,  "')", sep='')), 
              silent=TRUE)
    
    if(is(tr, "try-error")) { 
      if(verbose) cat("Adding", id, "failed -",
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


