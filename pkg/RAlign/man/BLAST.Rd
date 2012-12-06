\name{BLAST}
\alias{BLAST}
\title{
Basic Local Alignment Search Tool (BLAST)
}
\description{
Execute blastn from blast+.
}
\usage{
BLAST(x, db=NULL, BLAST_args="")
}
\arguments{
  \item{x}{ the query as an object of class XStringSet.} 
  \item{db}{ the database file to be searched (without file extension).}
  \item{BLAST_args}{ additional arguments in command-line style.}
}
\details{
    The blast+ software has to be installed and
    be executable by the system. Try \code{Sys.which("blastn")} to
    see if the program is properly installed.

    A database is needed. NCBI BLAST databases are updated daily 
    and may be downloaded via FTP 
    from \url{ftp://ftp.ncbi.nlm.nih.gov/blast/db/}.
}
\value{
Returns a data.frame containing the 
blast results.
}
\references{
BLAST+ \url{http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download}

}
\examples{
\dontrun{
seq <- readRNAStringSet(system.file("examples/RNA_example.fasta",
	package="RAlign"))

## shorten names
names(seq) <-  sapply(strsplit(names(seq), " "), "[", 1)
seq

## run BLAST (I have the NCBI database downloaded to ~/tmp/blast)
cl <- BLAST(seq, db="~/tmp/blast/refseq_rna")
cl[1:5, ]

}
}
\keyword{model}
