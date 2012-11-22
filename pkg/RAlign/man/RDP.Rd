\name{RDP}
\alias{RDP}
\alias{classifyRDP}
\alias{checkRDPOutput}
\alias{RDP_help}
\title{
Use RDP Classifier to classify DNA Sequences.
}
\description{
Uses the RDP classifier from qiime software to classify DNA sequences.
}
\usage{
classifyRDP(sequences=NULL, filepath=NULL)
checkRDPOutput(db,RDPfile="rdp/classify_tax_assignments.txt")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sequences}{an object of class DNAStringSet with the sequences to be classified.}
  \item{filepath}{Location of fasta file to be classified}
  \item{db}{GenDB containing details of sequences}
  \item{RDPfile}{Location of RDP classification file.}
}
\details{
The Qiime or Macqiime software has to be installed and 
be executable by the system. Try \code{Sys.which("macqiime")} or \code{Sys.which("qiime")} to 
see if the program is properly installed.
}
\value{
\code{classifyRDP} creates a classification file in the directory rdp.
\code{checkRDPOutput} checks the output of the RDP file against the data in the db.
}
\examples{
\dontrun{
### DNA
db<-openGenDB(system.file("examples/example.sqlite",package="RAlign"))
seq <- getSequences(db)
write.XStringSet(seq, file="seq.fasta")
classifyRDP(filepath="seq.fasta")
checkRDPOutput(db)
}
}
\keyword{model}
