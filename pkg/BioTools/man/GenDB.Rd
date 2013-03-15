\name{GenDB}
\alias{GenDB}
\alias{openGenDB}
\alias{createGenDB}
\alias{closeGenDB}
\alias{reopenGenDB}
\alias{listGenDB}
\alias{metaGenDB}
\alias{print.GenDB}
\title{Gene Sequence databases}
\description{
	Functions to create, open, close, and list tables from a 
	gene sequence database. 
}
\usage{
createGenDB(dbName, classification=GenClass16S_Greengenes(), drv=NULL, ...)
openGenDB(dbName, drv=NULL, ...)
reopenGenDB(db, ...)
closeGenDB(db)
listGenDB(db) 
metaGenDB(db, table=NULL) 
}
\arguments{
    \item{dbName}{Name of the database where the classification 
	and sequence information is to be stored.} 
    \item{db}{Handle of the database whose tables are to 
	be listed or which is to be closed.} 
    \item{drv}{Database driver, if none is provided then by 
	default SQLite driver is used.} 
    \item{classification}{Classification type object. By default, 
	a Greengenes 16S type classification is used. }
    \item{table}{Return meta data for a given table (NULL returns all meta data).}
    \item{...}{further arguments are passed on to the database driver.} 
}
\value{
    \code{createGenDB()}, \code{openGenDB()} and \code{reopenGenDB()} returns a handle to the database. 
    
    \code{listGenDB()} returns a vector of character strings.
}
\details{
	This group of functions are used to manage gene sequence database. 
	The function \code{createGenDB()} creates a
	database with three initial tables: sequences, classification, and
	metadata. The sequences table stores the sequences in character format
	with the organism name being the primary key. Similarly, the
	classification table stores the classification details such as kingdom,
	phylum, class, etc with the organism name being the primary key. The
	metadata table stores data about the tables in the database. 
	
	\code{openGenDB()} can be used to open an existing database using the database name (filename in SQLite),
	 \code{openGenDB()} can be used to open an existing database using the DB handle of a previousely closed GenDB,
	\code{listGenDB()} lists all the tables in an existing database,
	\code{metaGenDB()} lists the meta data about all the tables stored in
	the database, and \code{closeGenDB()} closes an open connection to a
	database.
}	
%\seealso{
%}
\examples{
db <- createGenDB("example.sqlite")
db

listGenDB(db)

### clean up
closeGenDB(db)
unlink("example.sqlite")
}
\keyword{database}
