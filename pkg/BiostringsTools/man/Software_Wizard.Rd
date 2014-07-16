\name{BiostringsTools_Software_Wizard}
\alias{BiostringsTools_Software_Wizard}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Software Installation Wizard for BiostringsTools 
}
\description{
This wizard helps with the installation of software used by BiostringsTools.
}
\usage{
BiostringsTools_Software_Wizard(RDP = FALSE, clustal = FALSE, 
  kalign = FALSE, MAFFT = FALSE, MUSCLE = FALSE, BLAST = FALSE, 
  BLAST16S = FALSE, boxshade = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{RDP, clustal, kalign, MAFFT, MUSCLE, BLAST, BLAST16S, boxshade}{ logical; 
  what would you like to install? Default is to install all except the BLAST16S
  data.
%%     ~~Describe \code{RDP} here~~
  }
}
\details{
For RDP you need a working Java installation. Note: Under Windows
you need to set PATH system variable so the java executible can be 
located.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\examples{
\dontrun{
BiostringsTools_Software_Wizard(RDP=TRUE)
}
}
%\keyword{ ~kwd1 }

