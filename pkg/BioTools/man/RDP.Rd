\name{RDP}
\alias{RDP}
\alias{predict}
\alias{predict.RDPClassifier}
\alias{print.RDPClassifier}
\alias{trainRDP}
\alias{removeRDP}
\alias{findAccuracy}
\title{
Ribosomal Database Project (RDP) Classifier for 16S rRNA
}
\description{
Use the RDP classifier to classify 16S rRNA sequences.
}
\usage{
RDP(dir = NULL)
\method{predict}{RDPClassifier}(object, newdata,
        confidence=.8, java_args="-Xmx1g", ...)
trainRDP(x, dir="classifier", rank="genus", java_args="-Xmx1g")
removeRDP(object)
findAccuracy(actual, predicted, rank)
}
\arguments{
  \item{dir}{directory where the classifier information is stored.} 
  \item{object}{a RDPClassifier object.}
  \item{newdata}{new data to be classified as a DNAStringSet.}
  \item{confidence}{numeric; minimum confidence level for classification. Results with lower confidence are replaced by NAs. Set to 0 to disable.}
  \item{java_args}{additional arguments for java (default sets the max. heap memory to 1GB).}
  \item{x}{an object of class DNAStringSet with the 16S rRNA sequences for training.}
  \item{actual}{A data frame containing the actual classification.}
  \item{predicted}{A data frame containing the predicted classification.}
  \item{rank}{Taxonomic rank at which the classification is learned/the accuracy is to be examined.}
  \item{...}{additional arguments (currently unused).}
}
\details{
Java and the RDP classifier has to be installed (See References for URL). 
You need to set an environment
variable called "RDP_JAR_PATH" with the path to the RDP classifier jar file.
You can set this variable using \code{Sys.setenv(RDP_JAR_PATH="/your_path/rdp_classifier-2.x.jar"} 
and test if it is set correctly with \code{Sys.getenv("RDP_JAR_PATH"}.
Note that on Windows the path will look like \code{"c:/your_path/rdp_classifier-2.x.jar"}.

\code{RDP()} creates a default classifier. Alternatively, a directory 
with the data for an existing classifier (created with \code{trainRDP()}) 
can be supplied.

\code{trainRDP()} creates a new classifier for the data in x and stores
the classifier information in \code{dir}. The data in
x needs to have annotations in the following format:

"<Kingdom>;<Phylum>;<Class>;<Order>;<Family>;<Genus">

A created classifier can be removed with \code{removeRDP()}. This will remove
the directory which stores the classifier information.
}
\value{
\code{RDP()} and \code{trainRDP()} return a RDPClassifier object.

\code{predict()} returns a data.frame containing the classification 
results for each sequence
(rows). The data.frame has an
attribure called "confidence" with a matrix containing the confidence values.
}
\references{
RDP Classifier \url{http://sourceforge.net/projects/rdp-classifier/}
}
\examples{
\dontrun{
### Use the default classifier
seq <- readRNAStringSet(system.file("examples/RNA_example.fasta",
	package="BioTools"))

## shorten names
names(seq) <-  sapply(strsplit(names(seq), " "), "[", 1)
seq

## use rdp for classification
predict(RDP(), seq)


### Train a custom RDP classifier on new data
trainingSequences <- readDNAStringSet(
    system.file("examples/trainingSequences.fasta", package="BioTools"))

customRDP <- trainRDP(trainingSequences)
customRDP

testSequences <- readDNAStringSet(
    system.file("examples/testSequences.fasta", package="BioTools"))
predict(customRDP, testSequences)

## clean up
removeRDP(customRDP)
}
}
\keyword{model}
