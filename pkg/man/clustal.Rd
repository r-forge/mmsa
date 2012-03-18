\name{clustal}
\alias{clustal}
\title{
Run Multiple Sequence Alignment (clustal) on a Set of Sequences
}
\description{
Executes clustal on a set of sequences to obtain a multiple sequence
alignment. 
}
\usage{
clustal(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
an object of class DNAStringSet with the sequences to be aligned.
}
}
\details{
The clustalw software has to be installed an be in the search path.
}
\value{
An object of class DNAMultipleAlignment (see \pkg{BioStrings}).
}
\references{
Larkin M., et al. Clustal W and Clustal X
version 2.0, Bioinformatics 2007 23(21):2947-29
}
\examples{
\dontrun{
### clustalW needs to be installed on your system!
### load 5 sequences
sequences <- read.DNAStringSet(
    system.file("examples/phylums/Proteobacteria.fasta",
            package="MMSA"),
    nrec=5)
sequences

### run clustalW
al <- clustal(sequences)
al
}
}
\keyword{manip}
