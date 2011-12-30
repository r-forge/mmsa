\name{addSequencesGreengenes}
\alias{addSequencesGreengenes}
\title{Reads a greengenes format fasta file}
\description{
	This function reads a fasta file from greengenes into the gene sequence database.
  }
\usage{
	addSequencesGreengenes(db, file)
}
\arguments{
  \item{db}{Name of an exisiting database where the classification and sequence information is to be stored} 
  \item{file}{Complete path of existing greengenes format fasta file containing classification and sequence information} 
}
\details{
	This function is used to read and parse an existing greengenes format fasta file and store it in a previously created database.
	As a precondition, SQLite should be installed and a database should already be created using the createGenDB function.
	The function reads the fasta file, parses the header information, and stores it in the classification table with the organism name
	as the Primary Key.
	The sequence information is stored in the sequence table with the organism name as the Primary Key.
	If an invalid database or file name is provided an error will be thrown.
	This function is used to open an existing SQLite database that is already present on the local or network database.

}
\seealso{
\code{\link{createGenDB}},
\code{\link{openGenDB}},
\code{\link{getSequences}}
}
\examples{
	db <- createGenDB("Firmicutes100.sqlite")
	addSequencesGreengenes(db, system.file("examples/Firmicutes100.fasta",package="MMSA"))
	db	
	closeGenDB(db)
	rm(db)
	unlink("Firmicutes100.sqlite")
}

