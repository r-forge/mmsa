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
s <- mutations(random_sequences(100), 100, .01, .01, .01)
s

dNSV <- NSVdist(s)
dNSV

### relationship with edit distance
dEdit <- stringDist(s)

df <- data.frame(dNSV=as.vector(dNSV), dEdit=as.vector(dEdit))
plot(sapply(df, jitter), cex=.1)
### add lower bound (2*k, for Manhattan distance)
abline(0,1/(2*3), col="red", lwd=2)
### add regression line
abline(lm(dEdit~dNSV, data=df), col="blue", lwd=2)
}
\keyword{model}
