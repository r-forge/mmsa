\name{dropNSVTable}
\alias{dropNSVTable}
\title{Deletes the specified NSV table from the database}
\description{
	Deletes the specified NSV format table from the database whose handle is provided as an argument.
  }
\usage{
	createNSVTable(db, tableName)
}
\arguments{
	\item{db}{Handle of the database} 
	\item{tableName}{Name of the NSV format table to be dropped}
}
\details{
	This function deletes the specified NSV format table from the database whose handle is provided as the first parameter.
	If a wrong handle or wrong table is provided, an error will be thrown.
}
\seealso{
\code{\link{createNSVTable}}
}
\references{

}
\examples{
SQLite()
drv<-dbDriver("SQLite");

## create/open GenDB
db <- createGenDB("16S.db")
addSequencesGreengenes(db, "Firmicutes100.fasta")
d<-getSequences(db, "Gen", "Bac")
createNSVTable(db, "NSV","kingdom","bacteria")
dropNSVTable(db,"NSV")
}
\keyword{manip}
