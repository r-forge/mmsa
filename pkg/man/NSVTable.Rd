\name{createNSVTable}
\alias{createNSVTable}
\alias{dropNSVTable}
\title{Manage NSV format tables}
\description{
	Create and Drop tables containing sequences in NSV format
  }
\usage{
	createNSVTable(db, tableName, rank=NULL, name=NULL, 
	    window=100, overlap=0, word=3, last_window=FALSE)
	dropNSVTable(db, tableName)
}
\arguments{
	\item{db}{Handle of the database} 
	\item{tableName}{Name of the new NSV table which will store the NSV format sequences.}
	\item{rank}{Name of Rank by which to filter. Fuzzy matching is performed if an incomplete name is specified.}
	\item{name}{Value of Rank by which to filter. Fuzzy matching is performed if an incomplete name is specified.}
	\item{window}{Value of window size.}
	\item{overlap}{Value of overlap parameter.}
	\item{word}{Value of word parameter.}
	\item{last_window}{Value of boolean last_window parameter.}
}
\details{
	\code{createNSVTable()} function converts character sequences into NSV format using the parameters specified or the default parameters.
	A new table is created in the database with the name equal to the second parameter.
	The new table is very similar to the sequences table in structure, but the sequence field contains NSVs instead of the character sequences.	
	Optional parameters whereRank and whereName parameters allow the conversion of only a subset of the entire sequences to NSV.

	\code{dropNSVTable()} function drops the specified table from the database whose handle is provided as the first argument.	
	If a wrong handle is provided, an error will be thrown.
}
\seealso{
\code{\link{getSequences}}
}
\examples{
	db<-createGenDB("MMSAExample.sqlite")
	addSequencesGreengenes(db, system.file("examples/phylums/Proteobacteria.fasta", package="MMSA"))
	
	### add the table with the counts to the db
	createNSVTable(db, "NSV")
	db

	### get NSVs for the first sequence
	getSequences(db, table="NSV")[1]

	### remove the table from the db
	dropNSVTable(db,"NSV")
	db
	closeGenDB(db)
	
	### clean up
	unlink("MMSAExample.sqlite")
		
}

