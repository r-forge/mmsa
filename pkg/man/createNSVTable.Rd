\name{createNSVTable}
\alias{createNSVTable}
\title{Creates a NSV Table by converting character sequences into NSV format}
\description{
	Created a NSV table by converting character sequences to NSV using the specified or default parameters.
  }
\usage{
	createNSVTable(db, tableName, whereRank, whereName, window=100, overlap=0, word=3, last_window=FALSE)
}
\arguments{
	\item{db}{Handle of the database} 
	\item{tableName}{Name of the new NSV table which will store the NSV format sequences.}
	\item{whereRank}{Name of Rank by which to filter. Fuzzy matching is performed if an incomplete name is specified. OPTIONAL}
	\item{whereName}{Value of Rank by which to filter. Fuzzy matching is performed if an incomplete name is specified. OPTIONAL}
	\item{window}{Value of window parameter, default value is 100. OPTIONAL}
	\item{overlap}{Value of overlap parameter, default value is 0. OPTIONAL}
	\item{word}{Value of word parameter, default value is 3. OPTIONAL}
	\item{last_window}{Value of boolean last_window parameter, default value is FALSE. OPTIONAL}
}
\details{
	This function converts character sequences into NSV format using the parameters specified or the default parameters.
	A new table is created in the database with the name equal to the second parameter.
	The new table is very similar to the sequences table in structure, but the sequence field contains NSVs instead of the character sequences.	
	Optional parameters whereRank and whereName parameters allow the conversion of only a subset of the entire sequences to NSV
	If a wrong handle is provided, an error will be thrown.
}
\seealso{
\code{\link{addSequencesGreengenes}},
\code{\link{nSequences}}
\code{\link{getSequences}}
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
}
\keyword{manip}
