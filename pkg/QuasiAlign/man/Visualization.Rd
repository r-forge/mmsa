\name{GenModel Visualization}
\alias{GenModelVisualization}
\alias{modelStatesPlot}
\alias{modelStatesBarPlot}
\alias{modelStatesLinePlot}
\alias{compareSequences}
\alias{findLargestCommon}
\title{
Visualize and compare sequences from EMM models.
}
\description{
Groups of functions to visualize, and compare sequences and NSVs from EMM models built with \code{saveClusterInfo=TRUE}.
}
\usage{
modelStatesPlot(model, states=NULL, numStates=5, ylab=TRUE, ...)
modelStatesBarPlot(model, ...)
modelStatesLinePlot(model, title=TRUE, ...)
compareSequences(model, sequences, ...)
findLargestCommon(model, limit=NULL)
}
\arguments{
	\item{model}{genModel object returned by the genModel and genModelDB classes.} 
	\item{states}{genModel states which are to be investigated.} 
	\item{numStates}{number of states to visualize.} 
	\item{sequences}{Gene sequences which are to be compared.}
	\item{limit}{Upper limit of number of sequences to be compared for common states.}
	\item{ylab}{Whether Y axis labels should be displayed.}
	\item{title}{Should the default title be shown for the modelStatesLine plot.}
	\item{...}{Other parameters to be passed to the plot function.}

}

\details{
	These groups of functions can be used to visually compare sequences and find common areas. These can provide a summary information of similar or conserved areas across sequences and are an alternative to Multiple Sequence Alignment which is computationally very expensive. 
}

\value{
\code{modelStatesPlot()}  returns a plot of the sequences highlighting
	which sequences and states are part of the given state(s). 

\code{modelStatesBarPlot()} returns a consensus-type bar plot showing which
state is dominant among the given window.

\code{modelStatesLinePlot()} returns a consensus-type line plot showing
    which state is dominant among the given window.

\code{compareSequences()} returns a plot of sequences highlighting the
	 given sequences and showing the common states and regions.

\code{findLargestCommon()} returns a list containing the most similar
     sequences from a model. Similar sequences are found by a count of the
     number of common states between them. The ith element of the list gives
     the i most similar sequences
	
}

\examples{
db <- openGenDB(system.file("examples/16S.sqlite",package="QuasiAlign"))

# model has to be built with saveClusterInfo=TRUE
model <- GenModelDB(db, rank="Gen", name="Syntro", table="NSV", 
    saveClusterInfo=TRUE) 
model

#create a model states plot for the 3 most preserved states
modelStatesPlot(model, numStates=3)

#visually compare sequences to find common areas
compareSequences(model, sequences=c(2,3,4,6))

#find the largest common sequences
findLargestCommon(model)

### clean up
closeGenDB(db)
}
%\keyword{model}
