#######################################################################
# BiostringsTools - Interfaces to several sequence alignment and 
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

### helper

.sq <-function(x) paste("'", x, "'", sep='')

.createdatatableGenDB <- function(db, table, type, annotation="") {
  if(length(grep(" ", table))) stop("Table name cannot contain spaces!")
  
  # data table stores the sequences as BLOB with id as PK
  seq <- "id TEXT PRIMARY KEY REFERENCES classification(id), data BLOB"
  tr <- try(
    dbSendQuery(db$db, 
                statement = paste("CREATE TABLE", .sq(table), 
                                  "(id TEXT PRIMARY KEY REFERENCES classification(id), data BLOB)")
                ), silent=TRUE)    
  if(is(tr, "try-error")) stop("Unable to create table ", table)
  
  #insert data into meta
  tr <- try(
    dbSendQuery(db$db,          
                statement = paste("INSERT INTO metaData  VALUES (", .sq(table),
                                  ",", .sq(type), ",", .sq(annotation),
                                  ")")
    ), silent=TRUE)
  if(is(tr, "try-error")) stop("Unable to write metaData for table", table)

}

createGenDB <- function(dbName, classification=GenClass16S_Greengenes(),
	 drv=NULL, ...) {

    if(file.exists(dbName)) 
	stop("GenDB already exist. Use openGenDB.\n")
    
    
    if(is.null(drv)) drv<-dbDriver("SQLite");
    db<-dbConnect(drv, dbname = dbName, ...);

    #first table stores the Classification with id as PK
    cl <- paste(.sq(names(classification)), "TEXT", collapse=', ')
    cl <- paste(cl, "PRIMARY KEY") ## the lowest rank is the primary key

    # classification table stores the class hierarchy
    tr <- try(
	    dbSendQuery(db, 
		    statement = paste("CREATE TABLE classification (",cl, ")"))
	    )
    if(is(tr, "try-error")) stop("Unable to create classification table")
    
   # metaData table stores the meta data
    meta<-"name TEXT, type TEXT, annotation TEXT"
    tr <- try(
      dbSendQuery(db, 
		    statement = paste("CREATE TABLE metaData (",
			    meta,  ")"))
	    )    
    if(is(tr, "try-error")) stop("Unable to create metaData table")
    
    ### start with no data table! 
    #.createdatatableGenDB(db, "sequences", "sequence")    

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

metaGenDB <- function(db, table=NULL) {
	meta <- dbGetQuery(db$db,"SELECT * from metaData")

	if(!is.null(table)) meta <- meta[meta[,"name"]==table, ]
	meta
    }

listGenDB <- function(db) setdiff(dbListTables(db$db), 
  c("classification", "metaData"))

dropTableGenDB <-  function(db, table) {

  if(any(tolower(table)==c("classification", "metaData"))) 
    stop("Cannot drop these tables!")
  dbSendQuery(db$db,
              statement = paste("DROP TABLE", .sq(table))
  )
  dbSendQuery(db$db,statement= paste("DELETE FROM metaData where name=",
                                     .sq(table),sep=''))
  dbCommit(db$db)
  invisible(NULL)
}

print.GenDB <- function(x, ...) {
    cat("Object of class GenDB with", nSequences(x, table="classification"), 
        "sequences\n")
    cat("DB File:", x$dbName, "\n")
    cat("Tables: ")
    cat(paste(listGenDB(x), collapse=", "), "\n")
}



