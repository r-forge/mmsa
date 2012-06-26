\name{Sequences}
\alias{Sequences}
\alias{getSequences}
\alias{addSequences}
\alias{nSequences}
\title{Genetic Sequences}
\description{
Genetic sequences are stored in a GenDB (as character strings with meta
information) or directly in R as objects of class DNAStringSet (package
\pkg{Biostrings}). \code{addSequences} adds sequences from a FASTA file to the
GenDB and \code{getSequences} retrieves sequences from the GenDB.
	
Sequences can also be represented in NSV format. 
See NSVSet for more information.
}
\usage{
	addSequences(db, file,
	    metaDataReader=GreengenesMetaDataReader, verbose=FALSE)
	getSequences(db, rank=NULL, name=NULL, table="sequences",
	    limit=NULL, random=FALSE, start=1, length=NULL,
	    partialMatch=TRUE)
	
	nSequences(db, rank=NULL, name=NULL)
}
\arguments{
    \item{db}{Handle of the GenDB database.} 
    \item{file}{Complete path of existing FASTA format file containing
	sequence information. }
    \item{metaDataReader}{function to process the meta data from the FASTA documentation line. }
    \item{verbose}{logical; print progress information.}
    \item{rank}{Name of Rank to be filtered by. For example: "kingdom", "phylum", etc.}
    \item{name}{Name of Rank to be filtered by. For example: if rank="kingdom", name can be "Bacteria".}
    \item{table}{Value of the NSV table which stores sequences in NSV format.}
    \item{limit}{Used to limit the number of sequnces that are considered for generating the model,by default all sequences are taken into account. Default value of -1 takes all sequences into account} 
    \item{random}{Get sequences in randomized order.} 
    \item{start}{Start position for the subsequence.} 
    \item{length}{Length of the subsequence (\code{NULL} means the whole sequence).} 
    \item{partialMatch}{logical; perform partial matching on the name.} 

}
\value{
	\code{getSequences()} returns an object containing the sequences.
	For table sequences the resulting object is a DNAStringSet for
	other tables it is a NSVSet object (see NSVSet).
}
\details{
	\code{addSequences()} reads and parses an existing FASTA format file
	and stores the sequences in an existing GenDB.  Meta data (biological
	classification) can also be extracted and stored in the GenDB for easy
	retrieval and selection of sequences. A custom reader for the meta data
	can be implemented and used.  The reader is called for each sequence,
	gets the content of the description line  and needs to return a vector
	containing the classification to be stored in the GenDB. The default is
	a reader for meta data used by the Greengenes project.
	
	\code{getSequences()} returns a list of sequences
	satisfying a search criteria specified by the optional rank and name
	parameters.  If no search criteria is specified, all the sequences are
	returned. By default, it searches the sequences table to obtain the
	list of sequences.  To obtain a list of NSVs, specify a NSV table as
	the "table" parameter. 
	
	\code{nSequences()} counts the total number of sequences
	present in a the sequence table of the database passed in as an
	argument.  If the optional rank and value parameters are supplied, the
	function filters the results based on the search criteria.
}
\seealso{
    \code{\link{createGenDB}},
    \code{\link{openGenDB}},
    \code{\link{GreengenesMetaDataReader}},
    \code{\link[Biostrings]{DNAStringSet}},
    \code{\link{NSVSet}}
}
\examples{
db <- createGenDB("example.sqlite")
addSequences(db,
    system.file("examples/phylums/Firmicutes.fasta",
        package="QuasiAlign"))


nSequences(db)

### select some sequences
s <- getSequences(db, rank="Class", name="Clostridia", limit=5)
s

### get the classification hierarchy for a the sequences
getHierarchy(db, "id", names(s))

### clean up
closeGenDB(db)
unlink("example.sqlite")
}
\keyword{database}
