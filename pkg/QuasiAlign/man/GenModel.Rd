\name{GenModel}
\alias{GenModel}
\alias{GenModelDB}
\alias{plot.GenModel}
\alias{print.GenModel}
\alias{prune,GenModel-method}
\alias{recluster}
\alias{getModelDetails}
\alias{getModelSequences}
\title{Genetic Sequences Models}
\description{
GenModel models sequences by clusters of similar NSVs and the transition
information between the clusters.
}
\usage{
GenModel(x, rank=NULL, name=NULL, measure="Manhattan", threshold=30, showClusterInfo=TRUE)

GenModelDB(db, rank=NULL, name=NULL, table ="NSV", measure="Manhattan", 
    threshold=30, selection=NULL, limit=NULL, random=FALSE, showClusterInfo=TRUE)

getModelDetails(model, state=NULL, db=NULL)
getModelSequences(db, model, state, table="sequences")

\S4method{prune}{GenModel}(x, ...)
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
	\item{random}{Used to select a random sample from the DB. Should be used for small selections.} 

	\item{...}{arguments passed on to the next method. For \code{prune}
	the arguments are passed on to \code{prune} in package \pkg{rEMM}. 
	The most important parameter is \code{count_threshold}. All
	clusters with a count of less or equal the threshold are pruned.}
}
\value{
	{\code{GenModel()}, \code{GenModelDB()} and \code{prune} return an object of class GenModel (a list  containing the model, rank, rankName, and number of sequences used in creating the model).}
	{\code{getModelDetails()} returns a list of states with the id of sequences in each state \code{getModelSequences()} returns a list of sequences in the model state specified by the modelState paramter.}

}
\details{
	genModel group of functions are  used after a gene sequence has been
	converted to NSV format.  After creating the model,
	\code{getModelDetails} and \code{getModelSequences} can be
	used to obtain meta information. 

	Also a plot function is provided for GenModel. Additional parameters
	are passed on to plot in package \pkg{rEMM}.
}
\seealso{\code{\link[rEMM]{prune}} in package \pkg{rEMM}}
\examples{
db<-createGenDB("example.sqlite")
processSequences(system.file("examples/phylums", package="QuasiAlign"), db)

### create a model for Genus Syntrophomonas
model <- GenModelDB(db, rank="Genus", name="Syntrophomonas") 
model

plot(model)


### look at sequence segments clustered in state 17
getModelDetails(model, 17)
s <- getModelSequences(db, model, state = 17)
s

nsv <- getModelSequences(db, model, state = 17, table="NSV")
nsv
plot(nsv)


### score 5 random sequence against the model
random_sequence <- getSequences(db, table="NSV", limit=5, random=TRUE)
scoreSequence(model, random_sequence)


### clean up
closeGenDB(db)
unlink("example.sqlite")
}
\keyword{model}
