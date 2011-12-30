\name{GenModel}
\alias{GenModel}
\alias{genModel}
\alias{plot.genModel}
\alias{plot}
\alias{classify}
\docType{data}
\title{
Manage Gene Sequence Models
}
\description{
Functions to generate, plot, and classify gene sequence models
}
\usage{
	genModel(db, rank, name, table, limit, measure, threshold, plus_one)
	plot.genModel(x, \dots)
}
\arguments{
	\item{db}{Handle of the database whose tables are to be listed or which is to be closed.} 
	\item{rank}{Name of Rank to be filtered by. For example: kingdom, phylum, etc. OPTIONAL}
	\item{name}{Name of Rank to be filtered by. For example: if rank=kingdom, name can be Bacteria. OPTIONAL}
	\item{table}{Name of the table containing NSV format sequences.} 
	\item{limit}{Used to limit the number of sequnces that are considered for generating the model,by default all sequences are taken into account. Default value of -1 takes all sequences into account} 
	\item{measure}{Character name of the dissimilarity measure used. See 'dist' in 'proxy' for available measures. Default is "Kullback"} 
	\item{threshold}{Numeric dissimilarity threshold used by the clustering algorithm for assigning a new observation to existing clusters. Default value is 0.10.} 
	\item{plus_one}{Boolean used to increase all the values in the NSV vectors by one. By default is TRUE.} 
	\item{x}{EMM model used in generic function for plotting.} 
	\item{...}{Further arguments passed on to 'plot.default'} 

}

\details{
	These group of functions are used after a gene sequence has been converted to NSV format. The first step is to create a gene sequence Model (genModel) using
	the genModel function, which returns the created emm. After that, this model can be used for plotting and classification purposes. 
}

\examples{
	db <- openGenDB(system.file("examples/16S.sqlite",package="MMSA"))
	createNSVTable(db, "NSV","kingdom","bacteria")
	emm<-genModel(db,rank="Gen",name="Thermo",table="NSV", measure="Kullback", threshold=0.10, plus_one=TRUE)
	plot.genModel(emm)
	dropNSVTable(db,"NSV")
	closeGenDB(db)
	rm(db)	
}
\keyword{datasets}
