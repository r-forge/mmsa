\name{GenModel}
\alias{GenModel}
\alias{GenModelDB}
\alias{plot.GenModel}
\alias{print.GenModel}
\alias{prune.GenModel}
\alias{recluster}
\alias{getClusteringDetails}
\alias{getClusteringSequences}
\title{
Create Models
}
\description{
Functions to create and handle models.
}
\usage{
GenModel(x, rank=NULL, name=NULL, measure="Manhattan", threshold=30)

GenModelDB(db, rank=NULL, name=NULL, table ="NSV", measure="Manhattan", 
    threshold=30, selection=NULL, limit=-1, showClusterInfo=TRUE)

getClusteringDetails(model, state=NULL)
getClusteringSequences(db, model, state, table="sequences")
}
\arguments{
	\item{x}{An object of class NSVSet.}
	\item{db}{Handle of the database whose tables are to be listed or which is to be closed.} 
	\item{rank}{Name of Rank to be filtered by. For example: "kingdom", "phylum", etc.}
	\item{name}{Name of Rank to be filtered by. For example: if rank="kingdom", name can be "Bacteria".}
	\item{table}{Name of the table containing NSV format sequences.} 
	\item{measure}{Character name of the dissimilarity measure used. See 'dist' in 'proxy' for available measures.} 
	\item{threshold}{Numeric dissimilarity threshold used by the clustering algorithm for assigning a new observation to existing clusters.}
	\item{selection}{Numeric vector indicating which sequences to take while creating the model. By default, all sequences are taken when constructing the model.} 
	\item{limit}{Used to limit the number of sequnces that are considered for generating the model,by default all sequences are taken into account. Default value of -1 takes all sequences into account.} 
	\item{model}{GenModel object returned by the GenModel and GenModelDB classes.} 
	\item{state}{state in the model about which details are required.} 
	\item{showClusterInfo}{Flag to indicate whether the model should contain cluster assignment information.} 

}
\value{
	{\code{GenModel()} and \code{GenModelDB()} return an object of class GenModel (a list  containing the model, rank, rankName, and number of sequences used in creating the model).}
	{\code{getClusteringDetails()} returns a list of states with the id of sequences in each state \code{getClusteringSequences()} returns a list of sequences in the model state specified by the modelState paramter.}

}
\details{
	genModel group of functions are  used after a gene sequence has been
	converted to NSV format.  After creating the model,
	\code{getClusteringDetails} and \code{getClusteringSequences} can be
	used to obtain meta information. 
}

\examples{
db<-createGenDB("example.sqlite")
processSequences(system.file("examples/phylums", package="MMSA"), db)

### create a model for Genus Syntrophomonas
model <- GenModelDB(db, rank="Genus", name="Syntrophomonas") 
model

### look at sequence segments clusterd in state 1
getClusteringDetails(model, 1)
getClusteringSequences(db, model, state = 1)

### score 5 random sequence against the model
random_sequence <- getSequences(db, table="NSV", limit=5, random=TRUE)
scoreSequence(model, random_sequence)

### visualize the model as a graph
plot(model)

### clean up
closeGenDB(db)
unlink("example.sqlite")
}
\keyword{model}
