\name{random_sequences}
\alias{random_sequences}
\title{
Create a Set of Random Sequences
}
\description{
Creates a set of random DNA, RNA or AA sequences.
}
\usage{
random_sequences(len, number=1, prob=NULL, type=c("DNA", "RNA", "AA"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{len}{sequence length}
  \item{number}{number of sequences in the set}
  \item{prob}{a vector with letter probabilities. 4 for DNA and RNA and 20
  for AA. The letter order is given by DNA_BASES, RNA_BASES and AA_ALPHABET 
  (first 20 letters). }
  \item{type}{sequence type}
}
%\details{
%}
\value{
A XStringSet.
}
\examples{
### create random sequences
seqs <- random_sequences(100, 10)
seqs

### calculate edit distance
stringDist(seqs)
}
\keyword{model}
