\name{16SrRNA}
\alias{16SrRNA}
\alias{16S}
\title{
RDP Ribosomal Database Project (RDP) 16S rRNA Training Set
}
\description{
RDP classifier data learned from the standard 16S rRNA training set.
}
%\usage{
%}
%\arguments{
%}
\details{
\code{rdp()} creates a default classifier using this data set.
}
%\value{
%}
\references{
RDP Classifier \url{http://sourceforge.net/projects/rdp-classifier/}
}
\examples{
### Use the default classifier
seq <- readRNAStringSet(system.file("examples/RNA_example.fasta",
	package="BiostringsTools"))

## shorten names
names(seq) <-  sapply(strsplit(names(seq), " "), "[", 1)
seq

## use rdp for classification 
pred <- predict(rdp(), seq)
pred
  
attr(pred, "confidence")  
}
\keyword{model}
