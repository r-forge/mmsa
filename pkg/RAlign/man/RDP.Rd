\name{RDP}
\alias{RDP}
\alias{classifyRDP}
\alias{RDP_help}
\title{
Ribosomal Database Project (RDP) Classifier for 16S rRNA
}
\description{
Use the RDP classifier to classify 16S rRNA sequences.
}
\usage{
classifyRDP(x, confidence=.8, property=NULL, java_args="-Xmx1g")
}
\arguments{
  \item{x}{an object of class DNAStringSet with the 16S rRNA sequences to be classified.}
  \item{confidence}{numeric; minimum confidence level for classification. Results with lower confidence are replaced by NAs. Set to 0 to disable.}
  \item{property}{a custom .property file.}
  \item{java_args}{additional arguments for java (default sets the max. heap memory to 1GB).}
}
\details{
Java and the RDP classifier has to be installed. You need to set an environment
variable called "RDP_JAR_PATH" with the path to the RDP classifier jar file.
You can set this variable using \code{Sys.setenv(RDP_JAR_PATH="/your_path/rdp_classifier-2.x.jar"} and test if it is set correctly with \code{Sys.getenv("RDP_JAR_PATH"}.
}
\value{
Returns a data.frame containing the classification results for each sequence
(rows). The data.frame has an
attribure called "confidence" with a matrix containing the confidence values.
}
\references{
RDP Classifier \url{http://sourceforge.net/projects/rdp-classifier/}
}
\examples{
\dontrun{
seq <- read.RNAStringSet(system.file("examples/RNA_example.fasta",
	package="RAlign"))

## shorten names
names(seq) <-  sapply(strsplit(names(seq), " "), "[", 1)
seq

## use rdp for classification
classifyRDP(seq)
}
}
\keyword{model}
