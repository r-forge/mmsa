\name{openGenDB}
\alias{openGenDB}
\title{Opens an existing SQLite database}
\description{
	This function opens an existing SQLite database whose name is provided as an argument.
  }
\usage{
db <- openGenDB("DBName.db")
}
\arguments{
  \item{DBName.db}{Name of the database to be opened} 
}
\details{
	This function is used to open an existing SQLite database that is already present on the local or network database.
	As a precondition, SQLite should be installed and running.
	The name of the database to be opened is given as the argument.
	If an invalid name is given, an error will be thrown. 
  	The database should be in standard MMSA format and should have tables for classification, sequences, and possibly one or more NSV tables.
	A meta table should store meta information about the table and the way they were created.
}
\seealso{
\code{\link{createGenDB}},
}
\references{

}
\examples{
SQLite()
drv<-dbDriver("SQLite");
db <- openGenDB("16S.db")
}
\keyword{manip}
