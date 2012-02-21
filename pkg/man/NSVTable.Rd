\name{createNSVTable}
\alias{createNSVTable}
\alias{dropNSVTable}
\title{Manage NSV format tables}
\description{
	Create and Drop tables containing sequences in NSV format
  }
\usage{
	createNSVTable(db, tableName, whereRank=NULL, whereName=NULL, window=100, overlap=0, word=3, last_window=FALSE)
	dropNSVTable(db, tableName)
}
\arguments{
	\item{db}{Handle of the database} 
	\item{tableName}{Name of the new NSV table which will store the NSV format sequences.}
	\item{whereRank}{Name of Rank by which to filter. Fuzzy matching is performed if an incomplete name is specified.}
	\item{whereName}{Value of Rank by which to filter. Fuzzy matching is performed if an incomplete name is specified.}
	\item{window}{Value of window size.}
	\item{overlap}{Value of overlap parameter.}
	\item{word}{Value of word parameter.}
	\item{last_window}{Value of boolean last_window parameter.}
}
\details{
	createNSVTable function converts character sequences into NSV format using the parameters specified or the default parameters.
	A new table is created in the database with the name equal to the second parameter.
	The new table is very similar to the sequences table in structure, but the sequence field contains NSVs instead of the character sequences.	
	Optional parameters whereRank and whereName parameters allow the conversion of only a subset of the entire sequences to NSV
	dropNSVTable function drops the specified table from the database whose handle is provided as the first argument.	
	If a wrong handle is provided, an error will be thrown.
}
\seealso{
\code{\link{addSequencesGreengenes}},
\code{\link{nSequences}},
\code{\link{getSequences}},
\code{\link{getRank}}
}
\examples{
	db<-createGenDB("MMSAExample.sqlite")
	addSequencesGreengenes(db,system.file("examples/phylums/Proteobacteria100.fasta",package="MMSA"))
	addSequencesGreengenes(db,system.file("examples/phylums/Planctomycetes100.fasta",package="MMSA"))
	
	createNSVTable(db, "NSV","kingdom","bacteria")
	
	dropNSVTable(db,"NSV")
	closeGenDB(db)
	unlink("MMSAExample.sqlite")
		
}

