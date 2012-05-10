\name{Hierarchy}
\alias{getRank}
\alias{getClassification}
\alias{getAllRanks}
\alias{getHierarchy}
\title{Classification Hierarchy for Sequences}
\description{
Inspect the classifiaction hierarchy of sequences in a GenDB.
Lists the distinct number of values at a given rank that are stored in the
database.  Also queries the classification hierarchy stored in the
classification table.  If no rank is provided, it defaults to the top-most
level i.e. kingdom.
  }
\usage{
getRank(db, rank, whereRank=NULL, whereName=NULL, all=FALSE)
getClassification(db)
getHierarchy(db, rank, name)
}
\arguments{
	\item{db}{Handle of the database} 
	\item{rank}{Name of the rank where the list of values is desired. Fuzzy matching is performed if an incomplete name is specified.}
	\item{whereRank}{Name of Rank by which to filter. Fuzzy matching is performed if an incomplete name is specified. }
	\item{whereName}{Value of Rank by which to filter. Fuzzy matching is performed if an incomplete name is specified. }
	\item{all}{logical; List the rank for each sequence separately?}
	\item{name}{Name of the rank at which the hierarchy object is required.}
}
\details{
	The \code{getRank()} function lists the distinct number of values at a
	given rank. For example, you might want to see how many distinct
	phylums are there in the table. As another example, if you would like
	to the distinct number of phylums where the kingdom is bacteria, you
	can use this function as:
	\code{getRank(db,"phylum","kingdom","bacteria")}.

	The \code{getClassification()} function lists the fields present in the
	classification table and can be used to find out the classification
	hierarchy.
	
	The \code{getHierarchy()} function returns a hierarchy object
	containing the hierarchy from the rank specified to the top of the
	classification tree.
}
%\seealso{
%}
\examples{
db <- openGenDB(system.file("examples/16S.sqlite",package="MMSA"))
db

### lists the distinct genuses
getRank(db, "genus")
table(getRank(db, "genus", all=TRUE))

### get the hierarchy object
getHierarchy(db, rank="genus", name="Syntrophomonas")

### get used classification hierarchy
getClassification(db)

closeGenDB(db)
}
\keyword{database}
