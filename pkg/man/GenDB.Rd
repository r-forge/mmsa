\name{GenDB}
\alias{openGenDB}
\alias{createGenDB}
\alias{closeGenDB}
\alias{listGenDB}
\alias{GenDB}
\alias{print.GenDB}
\docType{data}
\title{
	Manage Gene Sequence databases
}
\description{
	Functions to create, open, close, and list tables from a  gene sequence database on a local server. 
}
\usage{
	createGenDB(dbName, classification, type, drv)
	openGenDB(dbName, drv)
	closeGenDB(db)
	listGenDB(db) 
}
\arguments{
	\item{dbName}{Name of the database where the classification and sequence information is to be stored} 
	\item{db}{Handle of the database whose tables are to be listed or which is to be closed.} 
	\item{drv}{Database driver, if none is provided then by default SQLite driver is used. OPTIONAL.} 
	\item{classification}{Classification type object. By default, a Greengenes 16S type classification is used. OPTIONAL} 
	\item{type}{Type of annotation to be stored. By default, sequence annotations are stored. OPTIONAL} 

}
\value{
	\item{db}{createGenDB returns a handle to the database created}
}
\details{
	These group of functions are used to manage gene sequence database in a local database server. The function createGenDB
	creates a database with three initial tables: sequences, classification, and metadata. The sequences table stores the sequences in 
	character format with the organism name being the primary key. Similarly, the classification table stores the classification
	details such as kingdom, phylum, class, etc with the organism name being the primary key. The metadata table stores data about the
	tables in the database. 
	openGenDB can be used to open an existing database, listGenDB lists all the tables in an existing databases,  and closeGenDB closes an open connection to a database.
}
\seealso{
\code{\link{nSequences}},
\code{\link{getRank}},
\code{\link{addSequencesGreengenes}},
\code{\link{getSequences}},
\code{\link{getClassification}}
}
\examples{
	db <- createGenDB("GenDB.sqlite", drv=NULL, classification = GenClass16S_Greengenes(), type="sequence")
	openGenDB("GenDB.sqlite", drv=NULL)
	listGenDB(db)
	closeGenDB(db)
	rm(db)
	unlink("GenDB.sqlite")
	
}
\keyword{datasets}
