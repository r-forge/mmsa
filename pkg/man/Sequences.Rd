\name{Sequences}
\alias{Sequences}
\alias{getSequences}
\alias{nSequences}
\alias{count_sequences}
\alias{decodeSequence}
\title{
	Manage Gene sequences
}
\description{
	Functions to get, count, and decode sequences out of a gene sequence database.
}
\usage{
	nSequences(db, rank=NULL, name=NULL)
	getSequences(db, rank=NULL, name=NULL, table="sequences",
	    limit=-1, random=FALSE, start=1, length=100000000)
}
\arguments{
	\item{db}{Handle of the database} 
	\item{rank}{Name of Rank to be filtered by. For example: "kingdom", "phylum", etc.}
	\item{name}{Name of Rank to be filtered by. For example: if rank="kingdom", name can be "Bacteria".}
	\item{table}{Value of the NSV table which stores sequences in NSV format.}
	\item{limit}{Used to limit the number of sequnces that are considered for generating the model,by default all sequences are taken into account. Default value of -1 takes all sequences into account} 
	\item{random}{Get sequences in randomized order.} 
	\item{start}{Start position for the subsequence.} 
	\item{length}{Length of the subsequence.} 
	
}
\value{
	\code{getSequences()} returns a list of sequences.
}
\details{
	The function \code{getSequences()} returns a list of sequences satisfying a search criteria specified by the optional rank and name parameters.
	If no search criteria is specified, all the sequences are returned. By default, it searches the sequences table to obtain the list of sequences. 
	To obtain a list of NSVs, specify a NSV table as the "table" parameter. 
	
	The function \code{nSequences()} counts the total number of sequences present in a the sequence table of the database passed in as an argument.
	If the optional rank and value parameters are supplied, the function filters the results based on the search criteria.
}
\examples{
db <- openGenDB(system.file("examples/16S.sqlite",package="MMSA"))
db

nSequences(db,"class","Clostridia")
	
### get the first (limit=1) sequence
getSequences(db,"class","Clostridia",table="sequences", limit=1)
	
closeGenDB(db)
}
\keyword{manip}
