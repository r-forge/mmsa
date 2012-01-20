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
	genModel(db, rank, name, table, measure, threshold, plus_one, selection, limit)
}
\arguments{
	\item{db}{Handle of the database whose tables are to be listed or which is to be closed.} 
	\item{rank}{Name of Rank to be filtered by. For example: kingdom, phylum, etc. OPTIONAL}
	\item{name}{Name of Rank to be filtered by. For example: if rank=kingdom, name can be Bacteria. OPTIONAL}
	\item{table}{Name of the table containing NSV format sequences.} 
	\item{measure}{Character name of the dissimilarity measure used. See 'dist' in 'proxy' for available measures. Default is "Kullback"} 
	\item{threshold}{Numeric dissimilarity threshold used by the clustering algorithm for assigning a new observation to existing clusters. Default value is 0.10.} 
	\item{plus_one}{Boolean used to increase all the values in the NSV vectors by one. By default is TRUE.} 
	\item{selection}{Numeric vector indicating which sequences to take while creating the model. By default, all sequences are taken when constructing the model.} 
	\item{limit}{Used to limit the number of sequnces that are considered for generating the model,by default all sequences are taken into account. Default value of -1 takes all sequences into account} 

}
\value{
	\item{emm}{genModel function returns a list which has two parts - model and name. The model contains the actuall EMM model and the name contains annotations about the objects contained in the EMM.}
}
\details{
	These group of functions are used after a gene sequence has been converted to NSV format. The first step is to create a gene sequence Model (genModel) using
	the genModel function, which returns the created emm. After that, this model can be used for plotting and classification purposes. 
}

\examples{
	db <- openGenDB(system.file("examples/16S.sqlite",package="MMSA"))
	createNSVTable(db, "NSV","kingdom","bacteria")
	emm<-genModel(db,rank="Gen",name="Syntro",table="NSV", measure="Kullback", threshold=0.10, plus_one=TRUE)
	plot(emm$model, main=emm$name)
	dropNSVTable(db,"NSV")
	closeGenDB(db)
	rm(db)	
}
\keyword{datasets}
