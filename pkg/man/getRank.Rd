\name{getRank}
\alias{getRank}
\title{Lists the distinct values at a given rank}
\description{
	Lists the distinct number of values at a given rank that are stored in the database.
	If no rank is provided, it defaults to the top-most level i.e. kingdom.
  }
\usage{
getRank(db,rank,whereRank,whereName)
}
\arguments{
	\item{db}{Handle of the database} 
	\item{rank}{Name of the rank where the list of values is desired. Fuzzy matching is performed if an incomplete name is specified.}
	\item{whereRank}{Name of Rank by which to filter. Fuzzy matching is performed if an incomplete name is specified. OPTIONAL}
	\item{whereName}{Value of Rank by which to filter. Fuzzy matching is performed if an incomplete name is specified. OPTIONAL}
}
\details{
	This function lists the distinct number of values at a given rank. For example, you might want to see how many distinct phylums are
	there in the table. As another example, if you would like to the distinct number of phylums where the kingdom is bacteria, you can
	use this function as: getRank(db,"phylum","kingdom","bacteria")
	If the optional rank and value parameters are supplied, the function filters the results based on the search criteria.
	If a wrong handle is provided, an error will be thrown.
}
\seealso{
\code{\link{addSequencesGreengenes}},
\code{\link{getSequences}},
\code{\link{nSequences}},
\code{\link{getClassification}}
}
\examples{

	db <- openGenDB(system.file("examples/16S.sqlite",package="MMSA"))
	getRank(db,"genus")
	}
\keyword{manip}
