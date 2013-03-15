\name{clustal}
\alias{clustal}
\alias{clustal_help}
\title{
Run Multiple Sequence Alignment (clustal) on a Set of Sequences
}
\description{
Executes clustal on a set of sequences to obtain a multiple sequence
alignment. 
}
\usage{
clustal(x, param)
clustal_help()
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an object of class XStringSet (e.g., DNAStringSet) with the sequences to be aligned.}
  \item{param}{character string with the command line parameters for clustal (see output of \code{clustal_help()}).}

}
\details{
The clustalw software has to be installed and 
be executable by the system. Try \code{Sys.which("clustalw")} to 
see if the program is properly installed.
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
### DNA
dna <- readDNAStringSet(system.file("examples/DNA_example.fasta",
	package="BioTools"))
dna

al <- clustal(dna) 
al

### inspect alignment
detail(al)

### plot a sequence logo for the first 20 positions
plot(al, 1, 20)

### RNA
rna <- readRNAStringSet(system.file("examples/RNA_example.fasta",
	package="BioTools"))
rna 

al <- clustal(rna) 
al  

### Proteins 
aa <- readAAStringSet(system.file("examples/Protein_example.fasta",
	package="BioTools"))
aa 

al <- clustal(aa)
al
}
}
\keyword{model}
