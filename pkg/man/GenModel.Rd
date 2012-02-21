\name{GenModel}
\alias{GenModel}
\alias{genModel}
\alias{plot}
\docType{data}
\title{
Manage Gene Sequence Models
}
\description{
Functions to generate, plot, and classify gene sequence models
}
\usage{
	genModel(db, rank, name, table, measure="Kullback", threshold=0.10, plus_one=TRUE, selection=NULL, limit=-1)
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
	{genModel function returns a list  containing the model, rank, rankName, and number of sequences used in creating the model.}
}
\details{
	These group of functions are used after a gene sequence has been converted to NSV format. The first step is to create a gene sequence Model (genModel) using
	the genModel function, which returns the created emm. After that, this model can be used for plotting and classification purposes. 
}

\examples{
	db<-createGenDB("MMSAExample.sqlite")
	processSequencesGreengenes(system.file("examples/phylums",package="MMSA"),db)
    
	genModel(db,rank="Gen",name="Syntro",table="NSV", measure="Kullback", threshold=0.10, plus_one=TRUE)
	closeGenDB(db)
	unlink("MMSAExample.sqlite")
}
\keyword{datasets}
