\name{Sequences}
\alias{Sequences}
\alias{getSequences}
\alias{nSequences}
\alias{count_sequences}
\alias{decodeSequence}
\docType{data}
\title{
	Manage Gene sequences
}
\description{
	Functions to get, count, and decode sequences out of a gene sequence database.
}
\usage{
	nSequences(db,rank,name)
	getSequences(db,rank,name,table,limit,random)
}
\arguments{
	\item{db}{Handle of the database} 
	\item{rank}{Name of Rank to be filtered by. For example: kingdom, phylum, etc. OPTIONAL}
	\item{name}{Name of Rank to be filtered by. For example: if rank=kingdom, name can be Bacteria. OPTIONAL}
	\item{table}{Value of the NSV table which stores sequences in NSV format. If no value is provided, the "sequences" table is used as default.}
	\item{limit}{Used to limit the number of sequnces that are considered for generating the model,by default all sequences are taken into account. Default value of -1 takes all sequences into account} 
	\item{random}{Used to randomize the sequences before selection. Default value of -1 does not randomize sequences, a positive value will randomize sequences} 
	
}
\value{
	\item{d}{getSequences returns a list of sequences}
}
\details{
	The function getSequences returns a list of sequences satisfying a search criteria specified by the optional rank and name parameters.
	If no search criteria is specified, all the sequences are returned. By default, it searches the sequences table to obtain the list of sequences. 
	To obtain a list of NSVs, specify a NSV table as the "table" parameter. 
	The function nSequences counts the total number of sequences present in a the sequence table of the database passed in as an argument.
	If the optional rank and value parameters are supplied, the function filters the results based on the search criteria.
	If a wrong handle is provided, an error will be thrown.

}
\examples{
	db <- openGenDB(system.file("examples/16S.sqlite",package="MMSA"))
	nSequences(db,"class","Clostridia")
	getSequences(db,"class","Clostridia",table="sequences")
}
\keyword{datasets}
