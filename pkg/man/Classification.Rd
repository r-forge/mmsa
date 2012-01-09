\name{Classification}
\alias{Classification}
\alias{createModels}
\alias{processSequencesGreengenes}
\alias{classify}
\docType{data}
\title{
Generate models and classify sequences.
}
\description{
Function to generate and save EMM models from a given database and use them for classification of new sequences.
}
\usage{
processSequencesGreengenes(dir, db)
createModels(modelDir, rank, db)
classify(modelDir, seqFile)
}
\arguments{
	\item{db}{Handle of the database containing the sequences.}
	\item{dir}{Directory which contains the fasta format files which are to be read in.}
	\item{modelDir}{Directory where EMM models are stored.}
	\item{seqFile}{Name of the file containing sequences which are to be classified.}
	\item{rank}{Name of the rank (eg:phylum, class) at which the EMM model is to be built.}
}

\details{
	These group of functions are used to convert sequences into EMM models and then use these models for classification of new sequences.
	The function processSequencesGreengenes converts all the sequence files in a directory to NSV format in a database specified by the db parameter.
	The createModels function creates EMM models at a given rank (eg: phylum) and stores the models in a serialized format in the directory specified by the first parameter. 
	The function classify finds similarity between a sequence (whose classification hierarchy may be unknown) and the EMM models saved in a directory. The function goes through all the sequences in the seqFile (which should be in fasta format) and finds their similarity with all the models in the modelDir directory. The result is outputted in a tab separated file called ClassifyResults.txt.
}
\examples{
#create a db for demonstration purpose
db<-createGenDB("classificationExample.sqlite")
processSequencesGreengenes(system.file("examples/phylums",package="MMSA"),db)
createModels(system.file("examples/models",package="MMSA"),rank="phylum",db)
classify(system.file("examples/models",package="MMSA"),system.file("examples/testData.fasta",package="MMSA"))
#delete db to perform cleanup
unlink("classificationExample.sqlite")
}
\keyword{datasets}
