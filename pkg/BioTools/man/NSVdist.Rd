\name{NSVdist}
\alias{NSVdist}
\title{
Simple Distance between Sequences based on the k-mer Distribution 
}
\description{
Creates a k-mer frequency histogram for each sequence and then calculates
the difference between the histograms. The default setting creates
3-mer histograms and uses Manhattan distance, i.e., it is the 
number of 3-mers by which two sequences differ.
}
\usage{
NSVdist(x, k = 3, method = "Manhattan")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{  an object of class XStringSet containing the sequences.}
  \item{k}{ size of used k-mers.}
  \item{method}{ metric used to calculate the dissimilarity between two
  k-mer frequency distributions. }
}
%\details{
%}
\value{
A dist object.
}
%\references{
%}
\examples{
s <- mutations(random_sequences(100), 100)
s

dNSV <- NSVdist(s)
d

### relationship with edit distance
dEdit <- stringDist(s)

plot(jitter(dNSV/2/3), jitter(dEdit), xlim=c(0,10), ylim=c(0,10))
abline(0,1, col="red", lwd=2)
}
\keyword{model}
