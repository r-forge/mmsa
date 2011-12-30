\name{getClassification}
\alias{getClassification}
\title{Lists the classification for given database}
\description{
	Lists the classification hierarchy in a given database
  }
\usage{
getClassification(db)
}
\arguments{
	\item{db}{Handle of the database} 
}
\details{
	This function lists the fields present in the classification table of a given database provided by the argument.
	If a wrong handle is provided, an error will be thrown.
}
\seealso{
\code{\link{openGenDB}},
\code{\link{createGenDB}},
\code{\link{addSequencesGreengenes}}
}
\examples{
	db <- openGenDB(system.file("examples/16S.sqlite",package="MMSA"))
	getClassification(db)
}

