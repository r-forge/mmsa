\name{getSequences}
\alias{getSequences}
\title{Returns a list of Sequence satisfying a given search criteria}
\description{
	Returns a list of Sequences satisfying a given search criteria. If no search criteria is entered, the entire list of sequences is returned.
	By default, character sequences are returned. If you want NSV sequences, specify the name of the NSV table as the fourth argument.
  }
\usage{
d <- getSequences(db,whereRank,whereName) or
d <- getSequences(db,whereRank,whereName,tableName)
}
\arguments{
	\item{db}{Handle of the database} 
	\item{whereRank}{Name of Rank by which to filter. Fuzzy matching is performed if an incomplete name is specified. OPTIONAL}
	\item{whereName}{Value of Rank by which to filter. Fuzzy matching is performed if an incomplete name is specified. OPTIONAL}
	\item{tableName}{Value of the NSV table which stores sequences in NSV format. If no value is provided, the sequences table is used as default.}
}
\details{
	This function returns a list of sequences satisfying a search criteria specified by the optional second and third parameters.
	If no search criteria is specified, all the sequences are returned by default.
	By default, it searches the sequences table to obtain the list of sequences. To obtain a list of NSVs, specify a NSV table as the fourth parameter.
	If a wrong handle is provided, an error will be thrown.
}
\seealso{
\code{\link{getRank}},
\code{\link{nSequences}}
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

}
\keyword{manip}
