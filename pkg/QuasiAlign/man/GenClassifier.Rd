\name{GenClassifier}
\alias{GenClassifier}
\alias{print.GenClassifier}
\alias{predict.GenClassifier}
\alias{predict}
\alias{removeGenClassifier}
\alias{trainGenClassifier}
\docType{data}
\title{Generate Models to Classify Sequences}
\description{
Function to generate and save models which can be used to classify 
new sequences.
}
\usage{
GenClassifier(classifier_dir = NULL)
removeGenClassifier(object)
trainGenClassifier(classifier_dir, db, rank = "Phylum",
        table="NSV", selection=NULL, limit=NULL, add=FALSE, ...)
\method{predict}{GenClassifier}(object, newdata, rank="Phylum",
        method="supported_transitions", match_cluster="exact", ...)
}
\arguments{
	\item{classifier_dir}{directory to store the classifier data.}
	\item{object}{a GenClassifier object.}
	\item{rank}{Name of the rank (for example, phylum, class) at which the model is to be built.}
	\item{db}{Handle of the database containing the sequences.}
	\item{table}{a table with NSVs}
  \item{selection}{Numeric vector indicating which sequences to be used for creating the model. By default, all sequences are used for model creation.}
	\item{limit}{Number of sequences to consider while creating a genmodel. By default, all sequences are considered.}
	\item{add}{add new models to the classifer. Note: does not add sequences to existing 
  models bur replaces them.}
  \item{newdata}{Sequences to be classified as a NSVSet.}
	\item{method}{Method used to compute the similarity score between
	the model and a new sequence (see \code{score} in \pkg{rEMM} for available scoring methods).}
	\item{match_cluster}{How are NSVs matched against the model clusters ("exact", "nn" or "weighted").}
	\item{...}{additional arguments are passed on to \code{GenModelDB}.}
}
\value{\code{predict} returns a list with classification results.}
\seealso{
	\code{\link{GenModelDB}} and
	\code{\link[rEMM]{score}} in package \pkg{rEMM}.
}
\details{
\code{trainGenClassifier} produces GenModels 
at a specifier phylogenetic rank and stores them for later use.
}
\examples{
db <- openGenDB(system.file("examples/16S.sqlite",package="QuasiAlign"))

### create models at the phylum rank
modelRank <- "phylum"
classifier <- trainGenClassifier("models", db=db, rank=modelRank, limit=10)

### classify new sequences
sequences <- getSequences(db, modelRank, table="NSV", limit=5, random=TRUE)
sequences

pred <- predict(classifier, sequences, modelRank)
pred

### get true Phylum and create a confusion table
true <- getRank(db, rank="Phylum", whereRank="id", whereName=names(sequences))
true
table(pred$class, true)


### cleanup
closeGenDB(db)
removeGenClassifier(classifier)
}
\keyword{model}
