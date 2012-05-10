\name{NSVSet}
\alias{createNSVTable}
\alias{createNSVSet}
\alias{dropNSVTable}
\title{Create and Manage Sequences in NSV Format}
\description{
Functions to create, manage and store sequences in NSV (numerical sumarization vector) format. A NSV is a set of p-mer counts for a sequence cut into segments.

}
\usage{
	createNSVSet(x, 
	    window=100, overlap=0, word=3, last_window=FALSE)
	createNSVTable(db, tableName, rank=NULL, name=NULL, 
	    window=100, overlap=0, word=3, last_window=FALSE)
	dropNSVTable(db, tableName)
}
\arguments{
	\item{x}{An object of class DNAStringSet.} 
	\item{db}{Handle of the GenDB.} 
	\item{tableName}{Name of the new NSV table which will store the NSV format sequences.}
	\item{rank}{Name of Rank by which to filter. Fuzzy matching is performed if an incomplete name is specified.}
	\item{name}{Value of Rank by which to filter. Fuzzy matching is performed if an incomplete name is specified.}
	\item{window}{Value of window size.}
	\item{overlap}{Value of overlap parameter.}
	\item{word}{Value of word parameter.}
	\item{last_window}{Value of boolean last_window parameter.}
}
\details{
	\code{createNSVSet()} creates NSVs from the sequences stored in a
	DNAStringSet.
	
	\code{createNSVTable()} converts 
	sequences stored in the GenDB into NSV format 
	and stored the converted sequences in a new table in the GenDB.
	Optional parameters whereRank and whereName parameters allow 
	the conversion of only a subset of the entire sequences to NSV.

	\code{dropNSVTable()} function drops the specified table 
	from the database whose handle is provided as the first argument.	
}
%\seealso{
%}
\examples{
db<-createGenDB("example.sqlite")
addSequences(db, 
    system.file("examples/phylums/Proteobacteria.fasta", package="MMSA"))

### add the table with the counts to the db
createNSVTable(db, "NSV")
db

### get NSVs for the first sequence
getSequences(db, table="NSV")[1]

### remove the table from the db
dropNSVTable(db,"NSV")
db

### clean up
closeGenDB(db)
unlink("example.sqlite")
}
\keyword{database}
\keyword{model}

