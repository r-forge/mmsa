\name{listGenDB}
\alias{listGenDB}
\title{Lists the tables in a database}
\description{
	Lists the tables in a database. 
  }
\usage{
listGenDB(db)
}
\arguments{
	\item{db}{Handle of the database} 
}
\details{
	This function lists the number of tables in the database provided by the argument.
	If a wrong handle is provided, an error will be thrown.
}
\seealso{
\code{\link{openGenDB}},
\code{\link{createGenDB}},
}
\references{

}
\examples{
SQLite()
drv<-dbDriver("SQLite");

## create/open GenDB
db <- createGenDB("16S.db")
addSequencesGreengenes(db, "Firmicutes100.fasta")
listGenDB(db)

}
\keyword{manip}
