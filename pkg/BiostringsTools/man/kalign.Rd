\name{kalign}
\alias{kalign}
\alias{kalign_help}
\title{Multiple Sequence Alignment (kalign)}
\description{
Runs kalign progressive multiple sequence alignment on a set of sequences.
}
\usage{
kalign(x, param=NULL)
kalign_help()
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an object of class DNAStringSet with the sequences to be aligned.}
  \item{param}{character string with the command line parameters for kalign (see output of \code{kalign_help()}).}
}
\details{
The kalign software has to be installed (Debian/Ubuntu package: kalign) and 
be executable by the system. Try \code{Sys.which("kalign")} to 
see if the program is properly installed.
}
\value{
An object of class DNAMultipleAlignment (see \pkg{BioStrings}).
}
\references{
Lassmann T.,  Sonnhammer E. Kalign - an accurate and fast multiple 
sequence alignment algorithm, BMC Bioinformatics 2005, 6:298
}
\examples{
\dontrun{
BiostringsTools_Software_Wizard(kalign = TRUE)

dna <- readDNAStringSet(system.file("examples/DNA_example.fasta",
	package="BiostringsTools"))
dna

### align the sequences
al <- kalign(dna) 
al
}
}
