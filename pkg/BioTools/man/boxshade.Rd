\name{boxshade}
\alias{boxshade}
\alias{box}
\alias{boxshade_help}
\title{
Boxshade: Shading Multiple Aligned Sequences}
\description{
Executes boxshade on a multiple sequence alignment. 
}
\usage{
boxshade(x, file, dev="ps", param="-thr=0.5 -cons -ruler -def")
boxshade_help()
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a multiple alignment as an object of class 
  DNAMultipleAlignment, RNAMultipleAlignment or AAMultipleAlignment.}
  \item{file}{output file}
  \item{dev}{used output device. Available are: ps, eps, hpgl, rtf, crt, ansi, 
    vt, ascii, fig, pict and html.}
  \item{param}{character string with the command line parameters for clustal (see output of \code{boxshade_help()}).}

}
\details{
The boxshade software has to be installed (Debian/Ubuntu package: boxshade) and 
be executable by the system. Try \code{Sys.which("boxshade")} to 
see if the program is properly installed.
}
\value{
Only a file is created.
}
\references{
Boxshade has been written by Kay Hofmann and Michael D. Baron
}
\examples{
\dontrun{
rna <- readRNAStringSet(system.file("examples/RNA_example.fasta",
	package="BioTools"))
rna 

al <- clustal(rna) 

boxshade(al, file="alignment.ps")
}
}
\keyword{model}
