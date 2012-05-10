\name{addSequences}
\alias{addSequences}
\title{Add Sequences in FASTA Format to a GenDB}
\description{
	This function reads sequences in FASTA format,
	extracts meta information using a custom reader and
	adds the sequences to a genDB.
  }
\usage{
addSequences(db, file, 
    metaDataReader=GreengenesMetaDataReader, verbose=FALSE)
}
\arguments{
  \item{db}{Handle to an exisiting GenDB where the classification and 
  sequence information is to be stored.} 
  \item{file}{Complete path of existing FASTA format file containing 
  sequence information. } 
  \item{metaDataReader}{function to process the meta data from the FASTA documentation line. } 
  \item{verbose}{logical; print progress information.} 	
}
\details{
	This function is used to read and parse an existing FASTA format
	files and stores the sequences in an existing GenDB.  
	
	Meta data (biological classification) can also be extracted and stored
	in the GenDB for easy retrieval and selection of sequences. A custom
	reader for the meta data can be implemented and used. 
	The reader is called for each sequence, gets the content of the 
	description line  and needs to return a vector containing the 
	classification to be stored in the GenDB. The default is a reader 
	for meta data used by the Greengenes project.
}
\seealso{
\code{\link{createGenDB}},
\code{\link{openGenDB}},
\code{\link{GreengenesMetaDataReader}}
}
\examples{
db <- createGenDB("example.sqlite")
addSequences(db, 
    system.file("examples/phylums/Firmicutes.fasta",
    package="MMSA"))

db

### look at sequences
s <- getSequences(db)
s

### get the classification hierarchy for a sequence
### FIXME: this can be done better
getHierarchy(db, "Id",names(s[1]))
	
### clean up
closeGenDB(db)
unlink("example.sqlite")

### show the code of the meta data reader for Greengenes
GreengenesMetaDataReader
}
\keyword{database}
