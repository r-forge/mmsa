\name{GenModel}
\alias{GenModel}
\alias{genModel}
\alias{plot.GenModel}
\alias{print.GenModel}
\alias{prune.GenModel}
\alias{recluster.GenModel}
\alias{scoreSequence}
\docType{data}
\title{
Manage Gene Sequence Models
}
\description{
Functions to generate a model.
}
\usage{
genModel(db, rank=NULL, name=NULL, table, measure="Manhattan", 
    threshold=10, plus_one=TRUE, selection=NULL, limit=-1)
}
\arguments{
	\item{db}{Handle of the database whose tables are to be listed or which is to be closed.} 
	\item{rank}{Name of Rank to be filtered by. For example: "kingdom", "phylum", etc.}
	\item{name}{Name of Rank to be filtered by. For example: if rank="kingdom", name can be "Bacteria".}
	\item{table}{Name of the table containing NSV format sequences.} 
	\item{measure}{Character name of the dissimilarity measure used. See 'dist' in 'proxy' for available measures.} 
	\item{threshold}{Numeric dissimilarity threshold used by the clustering algorithm for assigning a new observation to existing clusters.}
	\item{plus_one}{Adds one to each transition count. This is equal to starting with a uniform prior for the transition count distribution, i.e. initially all transitions are equally likely. It also prevents the product of probabilities to be zero if a transition was never observed.}
	\item{selection}{Numeric vector indicating which sequences to take while creating the model. By default, all sequences are taken when constructing the model.} 
	\item{limit}{Used to limit the number of sequnces that are considered for generating the model,by default all sequences are taken into account. Default value of -1 takes all sequences into account.} 

}
\value{
	{\code{genModel()} function returns an object of class genModel (a list  containing the model, rank, rankName, and number of sequences used in creating the model).}
}
\details{
	This function is used after a gene sequence has been converted to NSV format. 
}

\examples{
	db<-createGenDB("example.sqlite")
	processSequences(system.file("examples/phylums", package="MMSA"), db)
	model <- genModel(db, rank="Gen", name="Syntro", table="NSV") 
	model
	#score sequence against the model
	random_sequence <- getSequences(db, table="NSV", limit=1, random=TRUE)[[1]]
	scoreSequence(model,random_sequence)

	plot(model)
	closeGenDB(db)
	
	### clean up
	unlink("example.sqlite")
}
\keyword{model}
