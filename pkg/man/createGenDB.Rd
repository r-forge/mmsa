\name{createGenDB}
\alias{createGenDB}
\title{Create a SQLite database}
\description{
  This function creates a SQLite database with the name provided and returns a handle to the database.
  }
\usage{
db <- createGenDB("DBName.db")
}
\arguments{
  \item{DBName.db}{Name of the database to be created} 
}
\details{
  	This function creates a database with the name supplied as the parameter on the local instance of SQLite.
	It creates three tables: one for storing classification information with the organism name as the Primary Key and 
	the ususal classification as the other columns eg: Kingdom, Phylum, Class,..
	The second table is an empty table which will be used to store the sequences read. It has the organism name as the
	Primary Key and will have another column to store the read sequences.
	In addition, a third meta table will be created which will store the meta information about the tables and the parameters
	used to create them. It returns a handle to the database created.
}
\seealso{
\code{\link{openGenDB}},
\code{\link{closeGenDB}}
}
\references{

}
\examples{
SQLite()
drv<-dbDriver("SQLite");

## create/open GenDB
db <- createGenDB("16S.db")
}
\keyword{manip}
