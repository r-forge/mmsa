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
    db <- list(
	    db=dbConnect(db$drv, dbname = db$dbName, ...), 
	    dbName=db$dbName, drv=db$drv)
    class(db) <- "GenDB"
    db
}


listGenDB <- function(db) setdiff(dbListTables(db$db), 
	c("classification", "metaData"))


metaGenDB <- function(db, table=NULL) {
	meta <- dbGetQuery(db$db,"SELECT * from metaData")

	if(!is.null(table)) meta <- meta[meta[,"name"]==table, ]
	meta
    }

print.GenDB <- function(x, ...) {
    cat("Object of class GenDB with", nSequences(x), "sequences\n")
    cat("DB File:", x$dbName, "\n")
    cat("Tables: ")
    cat(paste(listGenDB(x), collapse=", "), "\n")
}

