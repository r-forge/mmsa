\name{getRank}
\alias{getRank}
\alias{getClassification}
\alias{getHierarchy}
\title{Lists the distinct values at a given rank including the classification.}
\description{
	Lists the distinct number of values at a given rank that are stored in the database.
	Also queries the classification hierarchy stored in the classification table.
	If no rank is provided, it defaults to the top-most level i.e. kingdom.
  }
\usage{
getRank(db,rank,whereRank=NULL,whereName=NULL)
getHierarchy(db,rank,whereRank=NULL,whereName=NULL)
getClassification(db)
}
\arguments{
	\item{db}{Handle of the database} 
	\item{rank}{Name of the rank where the list of values is desired. Fuzzy matching is performed if an incomplete name is specified.}
	\item{whereRank}{Name of Rank by which to filter. Fuzzy matching is performed if an incomplete name is specified. }
	\item{whereName}{Value of Rank by which to filter. Fuzzy matching is performed if an incomplete name is specified. }
}
\details{
	The getRank function lists the distinct number of values at a given rank. For example, you might want to see how many distinct phylums are
	there in the table. As another example, if you would like to the distinct number of phylums where the kingdom is bacteria, you can
	use this function as: getRank(db,"phylum","kingdom","bacteria")
	The getHierarchy function is similar to getRank, except that it includes duplicates and lists all the values present for a given rank.
	If the optional rank and value parameters are supplied, the function filters the results based on the search criteria.
	Similarly, the whereRank and whereName parameters can be used to filter the results.
	The getClassification function lists the fields present in the classification table and can be used to find out the classification hierarchy.
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
	
	#lists the distinct genuses
	getRank(db,"genus")
	
	#lists all  the genuses
	getHierarchy(db,"genus")
	getClassification(db)
	}
\keyword{manip}
