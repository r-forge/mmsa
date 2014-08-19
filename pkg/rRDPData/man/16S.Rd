\name{16SrRNA}
\alias{16SrRNA}
\alias{16S}
\title{
RDP Classifier Model for RDP 16S rRNA Trainingset No 9
}
\description{
Model data for the RDP classifier trained with the RDP 16S rRNA 
Trainingset No 9.
This is the default model data used by RDP.
}
%\usage{
%}
%\arguments{
%}
\details{
\code{rdp()} creates a default classifier using this data set.

The model data 
is stored in \code{system.file("16srrna", package="rRDPData")} and
consists of several files in proprietary format used by the RDP software. 
Since RDP is a naive Bayes classifier the data mainly contain conditional 
word (8-mer) probabilities.
}
%\value{
%}
\references{
RDP Classifier \url{http://sourceforge.net/projects/rdp-classifier/}
}
\examples{
### Use the RDP with the default classifier model data
seq <- readRNAStringSet(system.file("examples/RNA_example.fasta",
	package="rRDP"))

## shorten names
names(seq) <-  sapply(strsplit(names(seq), " "), "[", 1)
seq

## use RDP for classification 
pred <- predict(rdp(), seq)
pred
  
attr(pred, "confidence")  
}
\keyword{model}
