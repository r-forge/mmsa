\name{Visualization}
\alias{modelStatesPlot}
\alias{compareSequences}
\alias{findLargestCommon}
\title{
Visualize and compare sequences from EMM models.
}
\description{
Groups of functions to visualize, and compare sequences and NSVs from EMM models.
}
\usage{
modelStatesPlot(model, states)
compareSequences(model, sequences)
findLargestCommon(model, limit=NULL)
}
\arguments{
	\item{model}{genModel object returned by the genModel and genModelDB classes.} 
	\item{states}{genModel states which are to be investigated.} 
	\item{sequences}{Gene sequences which are to be compared.}
	\item{limit}{Upper limit of number of sequences to be compared for common states.}

}
\value{
	{\code{modelStatesPlot()}  returns a plot of the sequences highlighting which sequences and states are part of the given state(s). 
     \code{compareSequences()} returns a plot of sequences highlighting the given sequences and showing the common states and regions
     \code{findLargestCommon()} returns a list containing the most similar sequences from a model. Similar sequences are found by a count of the number
	of common states between them. The ith element of the list gives the i most similar sequences}

}

\details{
	These groups of functions can be used to visually compare sequences and find common areas. These can provide a summary information of similar or conserved areas across sequences and are an alternative to Multiple Sequence Alignment which is computationally very expensive. 
}

\examples{
db<-createGenDB("example.sqlite")
processSequences(system.file("examples/phylums", package="MMSA"), db)
model <- GenModelDB(db, rank="Gen", name="Syntro", table="NSV") 
model

#create a model states plot for multiple states
modelStatesPlot(model,c(17,18,19))
#visually compare sequences to find common areas
compareSequences(model,c(2,3,4,6))
#find the largest common sequences
common <- findLargestCommon(model)


### clean up
closeGenDB(db)
unlink("example.sqlite")
}
%\keyword{model}
