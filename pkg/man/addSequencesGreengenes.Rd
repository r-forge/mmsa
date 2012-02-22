\name{addSequencesGreengenes}
\alias{addSequencesGreengenes}
\title{Reads a greengenes format FASTA file}
\description{
	This function reads a FASTA file from greengenes into the gene sequence database.
  }
\usage{
	addSequencesGreengenes(db, file, verbose=FALSE)
}
\arguments{
  \item{db}{Handle to an exisiting database where the classification and sequence information is to be stored} 
  \item{file}{Complete path of existing greengenes format fasta file containing classification and sequence information} 
  \item{verbose}{If true, the function will print out a statement after reading every 100 sequences.} 	
}
\details{
	This function is used to read and parse an existing greengenes format FASTA file and store it in a previously created database.
	As a precondition, SQLite should be installed and a database should already be created using the createGenDB function.
	The function reads the FASTA file, parses the header information, and stores it in the classification table with the organism name
	as the Primary Key.
	The sequence information is stored in the sequence table with the organism name as the Primary Key.
	If an invalid database or file name is provided an error will be thrown.


}
\seealso{
\code{\link{createGenDB}},
\code{\link{openGenDB}},
\code{\link{getSequences}}
}
\examples{
db <- createGenDB("Firmicutes.sqlite")
addSequencesGreengenes(db, system.file("examples/phylums/Firmicutes.fasta",package="MMSA"))

db

closeGenDB(db)
	
### clean up
unlink("Firmicutes.sqlite")
}

