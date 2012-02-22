\name{GenDB}
\alias{openGenDB}
\alias{createGenDB}
\alias{closeGenDB}
\alias{listGenDB}
\alias{GenDB}
\alias{print.GenDB}
\title{
	Manage Gene Sequence databases
}
\description{
	Functions to create, open, close, and list tables from a  gene sequence database on a local server. 
}
\usage{
	createGenDB(dbName, classification=GenClass16S_Greengenes(), drv=NULL)
	openGenDB(dbName, drv=NULL)
	closeGenDB(db)
	listGenDB(db) 
}
\arguments{
	\item{dbName}{Name of the database where the classification and sequence information is to be stored} 
	\item{db}{Handle of the database whose tables are to be listed or which is to be closed.} 
	\item{drv}{Database driver, if none is provided then by default SQLite driver is used.} 
	\item{classification}{Classification type object. By default, a Greengenes 16S type classification is used. } 

}
\value{
	{\code{createGenDB()} and \code{openGenDB()} returns a handle to the database. \code{listGenDB()} returns a vector of character strings.}
}
\details{
	These group of functions are used to manage gene sequence database in a local database server. The function \code{createGenDB()}
	creates a database with three initial tables: sequences, classification, and metadata. The sequences table stores the sequences in 
	character format with the organism name being the primary key. Similarly, the classification table stores the classification
	details such as kingdom, phylum, class, etc with the organism name being the primary key. The metadata table stores data about the
	tables in the database. 
	
	\code{openGenDB()} can be used to open an existing database, \code{listGenDB()} lists all the tables in an existing databases,  and \code{closeGenDB()} closes an open connection to a database.
}
%\seealso{
%}
\examples{
db <- createGenDB("example.sqlite")
db

listGenDB(db)

closeGenDB(db)
	
### clean up
unlink("example.sqlite")
	
}
\keyword{database}
