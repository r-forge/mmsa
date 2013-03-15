\name{simRank}
\alias{simRank}
\title{
Compute the SimRank Similarity between Sets of Sequences
}
\description{
Computes the SimRank similarity (number of shared unique k-mers over the 
smallest number of unique k-mers.)
}
\usage{
simRank(x, k = 7)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{  an object of class DNAStringSet containing the sequences.
}
  \item{k}{ size of used k-mers.
}
}
%\details{
%}
\value{
A similarity object of class "simil" (see \pkg{proxy}).
}
\references{
Santis et al, Simrank: Rapid and sensitive general-purpose k-mer 
search tool, BMC Ecology 2011, 11:11
}
\examples{
### load sequences
sequences <- readDNAStringSet(
    system.file("examples/DNA_example.fasta", 
	package="BioTools"))
sequences

### compute similarity
simil <- simRank(sequences)

### use hierarchical clustering
hc <- hclust(pr_simil2dist(simil))
plot(hc)
}
\keyword{model}
