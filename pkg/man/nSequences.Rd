\name{nSequences}
\alias{nSequences}
\title{Counts the number of sequences in a database}
\description{
	Counts the total number of sequences which have been read into a database.
	Can provide Optional parameters for Rank and Value which will limit the number of results to those that match the criterion.
  }
\usage{
nSequences(db)
}
\arguments{
	\item{db}{Handle of the database} 
	\item{rank}{Name of Rank to be filtered by. For example: kingdom, phylum, etc. OPTIONAL}
	\item{value}{Value of Rank to be filtered by. For example: if rank=kingdom, value can be Bacteria. OPTIONAL}
}
\details{
	This function counts the total number of sequences present in a the sequence table of the database passed in as an argument.
	If the optional rank and value parameters are supplied, the function filters the results based on the search criteria.
	If a wrong handle is provided, an error will be thrown.
}
\seealso{
\code{\link{openGenDB}},
\code{\link{getSequences}},
}
\references{

}
\examples{
SQLite()
drv<-dbDriver("SQLite");

## create/open GenDB
db <- createGenDB("16S.db")
addSequencesGreengenes(db, "Firmicutes100.fasta")
nSequences(db)

}
\keyword{manip}
