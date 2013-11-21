\name{GenDB_Sequences}
\alias{GenDB_Sequences}
\alias{Sequences}
\alias{getSequences}
\alias{addSequences}
\alias{nSequences}
\alias{Annotation_Id}
\alias{Annotation_Greengenes}
\alias{Annotation_RDP}
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
addSequences(db, sequences, table="sequences",
	    annotation=Annotation_Greengenes, verbose=FALSE)
	
getSequences(db, rank=NULL, name=NULL, table="sequences",
	    limit=NULL, random=FALSE, start=NULL, length=NULL,
	    removeUnknownSpecies=FALSE, annotation=Annotation_Id)

nSequences(db, rank=NULL, name=NULL, table="classification")

Annotation_Id(annotation, decode)
Annotation_Greengenes(annotation, decode=TRUE)  
Annotation_RDP(annotation, decode=TRUE)
}
\arguments{
    \item{db}{Handle of the GenDB database.} 
    \item{sequences}{Complete path of existing FASTA format file containing
	sequence information, a directory name with FASTA format files, or a DNAStringSet
    with sufficient annotations in the name attribute.}
    \item{verbose}{logical; print progress information.}
    \item{rank}{Name of Rank to be filtered by. For example: "phylum", etc.}
    \item{name}{Name of Rank to be filtered by. For example: if rank="kingdom", name can be "Bacteria".}
    \item{table}{Value of the table which stores sequences.}
    \item{limit}{Used to limit the number of sequences that are considered for generating the model,by default all sequences are taken into account. Default value of -1 takes all sequences into account} 
    \item{random}{Get sequences in randomized order.} 
    \item{start}{Start position for the subsequence.} 
    \item{length}{Length of the subsequence (\code{NULL} means the whole sequence).} 
    \item{removeUnknownSpecies}{logical; whether unknown species should be removed from the result.} 
    \item{annotation}{Annotation format. Options are Annotation_Id, 
                      Annotation_Greengenes and Annotation_RDP. For the Annotation_x
                      functions, annotation contains the FASTA file description 
                      (obtained from a DNAStringSet via name).}
          
    \item{decode}{Decode or encode an annotation?}
}
\value{
	\code{getSequences()} returns an object containing the sequences
	as a DNAStingSet.
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
	returned. Partial matching is used to match the rank name to the known 
   rank names (e.g., P is matched to Phylum). Partial matching for name
   is available with the wild card character % (e.g., De% matches all names
   starting with De). 
	
	\code{nSequences()} counts the total number of sequences
	present in a the sequence table of the database passed in as an
	argument.  If the optional rank and value parameters are supplied, the
	function filters the results based on the search criteria.
}
\seealso{
  \code{\link{GenDB}},
  \code{\link{GenDB_Hierarchy}},
  \code{\link{GenDB_Greengenes}},
  \code{\link[Biostrings]{DNAStringSet}}
}
\examples{
db <- createGenDB("example.sqlite")
addSequences(db,
    system.file("examples/Firmicutes.fasta",
        package="BioTools"))

metaGenDB(db)
nSequences(db)

### select some sequences
s <- getSequences(db, rank="Class", name="Clostridia", limit=5)
s

### get the classification hierarchy for a the sequences
getHierarchy(db, "id", names(s))

### get sequences with other annotation formats
s <- getSequences(db, rank="Class", name="Clostridia", limit=5, 
                  annotation=Annotation_Greengenes)
s
names(s)

s <- getSequences(db, rank="Class", name="Clostridia", limit=5,
                  annotation=Annotation_RDP)
s
names(s)

### partial matching (G matches Genus and De% matches all names starting with De)
s <- getSequences(db, rank="G", name="De\%")
s
getRank(db, rank="Genus", whereRank="id", whereName=names(s), count=TRUE)

### clean up
closeGenDB(db)
unlink("example.sqlite")
}
\keyword{database}
