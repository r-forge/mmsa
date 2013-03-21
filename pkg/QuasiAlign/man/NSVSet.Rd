\name{NSVSet}
\alias{NSVSet}
\alias{createNSVTable}
\alias{createNSVSet}
\alias{dropNSVTable}
\alias{getNSVs}
\alias{getSequences}
\alias{plot.NSVSet}
\alias{c.NSVSet}
\alias{[.NSVSet}
\alias{processSequences}
\title{Sequences in NSV Format}
\description{
Sequences 
can be represented in NSV (numerical sumarization vector) format. 
A NSV is a set of p-mer counts for a sequence cut into segments.

NSVSets can be created from DNAStringSet with \code{createNSVSet}
or directly from data stored in a GenDB with \code{createNSVTable}.
}
\usage{
createNSVSet(x, window=100, overlap=0, word=3, 
    last_window=FALSE, allOffsets=FALSE)

createNSVTable(db, table="NSV", rank=NULL, name=NULL, 
    window=100, overlap=0, word=3, last_window=FALSE, 
    limit=NULL, removeUnknownSpecies=FALSE, allOffsets=FALSE)

dropNSVTable(db, table)

getNSVs(db,  rank=NULL, name=NULL,
    table="NSV", limit=NULL, random=FALSE, start=1, length=NULL,
    partialMatch=TRUE, removeUnknownSpecies=FALSE, annotation="id")

getSequences(db,  rank=NULL, name=NULL,
    table="sequences", limit=NULL, random=FALSE, start=1, length=NULL,
    partialMatch=TRUE, removeUnknownSpecies=FALSE, annotation="id")

processSequences(dir, db, metaDataReader=GreengenesMetaDataReader, ...)
}
\arguments{
	\item{x}{An object of class DNAStringSet.} 
	\item{db}{Handle of the GenDB.} 
	\item{table}{Name of the new NSV table which will store the NSV format sequences.}
	\item{rank}{Name of Rank by which to filter. Fuzzy matching is performed if an incomplete name is specified.}
	\item{name}{Value of Rank by which to filter. Fuzzy matching is performed if an incomplete name is specified.}
	\item{window}{Value of window size.}
	\item{overlap}{Value of overlap parameter.}
	\item{word}{Value of word parameter.}
	\item{last_window}{Value of boolean last_window parameter.}
	\item{limit}{Length of sequences to limit.}
	\item{removeUnknownSpecies}{Logical. Should unknown species be not be included while creating NSV Table?}
	\item{allOffsets}{Logical. Should NSV be created using all possible offsets?}
	\item{random}{Get sequences in randomized order.}
	\item{start}{Start position for the subsequence.}
	\item{length}{ignored.}
	\item{partialMatch}{logical; perform partial matching on the name for rank.}
	\item{annotation}{Annotation of the returned sequences. Options are Greengenes id or RDP format.}
	\item{dir}{Directory which contains the FASTA format files which are to be read in.}
	\item{metaDataReader}{Reader function used to interpret the FASTA meta information (see \pkg{BioTools}). }
	\item{...}{additional aruments are passed on. }
}
\details{
	\code{createNSVSet()} creates NSVs from the sequences represented as a
	DNAStringSet.
	
	\code{createNSVTable()} converts 
	sequences stored in the GenDB into NSV format 
	and stored the converted sequences in a new table in the GenDB.
	Optional parameters whereRank and whereName parameters allow 
	the conversion of only a subset of the entire sequences to NSV.

	\code{dropNSVTable()} function drops the specified table 
	from the database whose handle is provided as the first argument.	

	NSVSets stored in a GenDB can be retrieved using \code{getNSVs} 
	where table is
	the name of the table used in \code{createNSVTable}.
	For convenience we also provide a wrapper for 
	\code{getSequences} which will return sequences or NSVs depending
	on the given table.

	\code{plot} creates a barplot with the average p-mer counts
	in the NSVSet.

	\code{processSequences} is a convenience function that 
	loads all the sequence files in a
	directory into the GenDB and converts them into
	NSV format.
}
%\seealso{
%   \code{\link{getSequences}} 
%}
\examples{
db<-createGenDB("example.sqlite")
addSequences(db, 
    system.file("examples/phylums/Proteobacteria.fasta", 
	package="QuasiAlign"))

## add the table with the counts for all sequences in the GenDB
## we convert ony 5 sequences (limit) to save time
createNSVTable(db, table="NSV", limit=10)
db

### get NSVs for the first sequence
s <- getSequences(db, table="NSV", limit=5)
s

### select the first element (sequence)
s[1]

### select the first segment in all sequences
s[, 1]

### plot the p-mer distribution for the first segments
### (average counts with whiskers show the min/max counts)
plot(s[, 1])

### remove the table from the db
dropNSVTable(db,"NSV")
db

### clean up
closeGenDB(db)
unlink("example.sqlite")
}
\keyword{database}
\keyword{model}

