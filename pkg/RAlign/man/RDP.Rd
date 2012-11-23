\name{RDP}
\alias{RDP}
\alias{classifyRDP}
%\alias{checkRDPOutput}
\alias{RDP_help}
\title{
Use RDP Classifier to classify DNA Sequences.
}
\description{
Uses the RDP classifier from qiime software to classify DNA sequences.
}
\usage{
classifyRDP(x,param=NULL)
}
%checkRDPOutput(db,RDPfile="rdp/classify_tax_assignments.txt")
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an object of class DNAStringSet with the sequences to be classified.}
  \item{param}{a list of additional parameters (currently unused).}
%  \item{filepath}{Location of fasta file to be classified}
%  \item{db}{GenDB containing details of sequences}
%  \item{RDPfile}{Location of RDP classification file.}
}
\details{
The Qiime or Macqiime software has to be installed and 
be executable by the system. Try \code{qiime print_qiime_config} or \code{macqiime print_qiime_config} to 
see if the program is properly installed.
Note: Qiime is currently not available for Windows.
}
\value{
Returns a data.frame containing the classification.
%\code{classifyRDP} creates a classification file in the directory rdp.
%\code{checkRDPOutput} checks the output of the RDP file against the data in the db.
}
\references{
RDP Classifier \code{http://sourceforge.net/projects/rdp-classifier/},
Qiime \code{http://qiime.org/}
}

\examples{
\dontrun{
seq <- read.RNAStringSet(system.file("examples/RNA_example.fasta",
	package="RAlign"))

## shorten names
names(seq) <-  sapply(strsplit(names(seq), " "), "[", 1)
seq

## use rdp for classification
classifyRDP(seq)
}
}
\keyword{model}
