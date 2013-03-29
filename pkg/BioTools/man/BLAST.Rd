\name{BLAST}
\alias{BLAST}
\alias{BLAST_help}
\alias{predict.BLAST}
\alias{print.BLAST}
\title{
Basic Local Alignment Search Tool (BLAST)
}
\description{
Execute blastn from blast+.
}
\usage{
BLAST(db=NULL)
BLAST_help()
\method{predict}{BLAST}(object, newdata,
        BLAST_args="",  ...)
\method{print}{BLAST}(x, info=FALSE, ...)
}
\arguments{
  \item{db}{ the database file to be searched (without file extension).}
  \item{object, x}{ A BLAST object. }
  \item{newdata}{ the query as an object of class XStringSet.} 
  \item{BLAST_args}{ additional arguments in command-line style.}
  \item{info}{ print information about the database (needs blastdbcmd 
  in the path).}
  \item{...}{ additional arguments are ignorred.}
}
\details{
    The blast+ software has to be installed and
    be executable by the system. Try \code{Sys.which("blastn")} to
    see if the program is properly installed (
    Debian/Ubuntu package: ncbi-blast+).

    A database is needed. NCBI BLAST databases are updated daily 
    and may be downloaded via FTP 
    from \url{ftp://ftp.ncbi.nlm.nih.gov/blast/db/}.
}
\value{
\code{BLAST} returns a blast database object which can be used for queries 
(via \code{predict}). \code{predict}
returns a data.frame containing the 
blast results.
}
\references{
BLAST+ \url{http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download}

}
\examples{
\dontrun{
seq <- readRNAStringSet(system.file("examples/RNA_example.fasta",
	package="BioTools"))

## shorten names
names(seq) <-  sapply(strsplit(names(seq), " "), "[", 1)
seq

## load a BLAST database (replace db with the location + name of the BLAST DB)
blast <- BLAST(db="~/tmp/blast/16SMicrobial")
blast

print(blast, info=TRUE)

## query a sequence using BLAST
cl <- predict(blast, seq[1,])
cl[1:5,]

}
}
\keyword{model}
