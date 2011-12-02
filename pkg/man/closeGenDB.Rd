\name{closeGenDB}
\alias{closeGenDB}
\title{Closes an open SQLite database}
\description{
	This function closes an existing open SQLite database. 
  }
\usage{
db <- closeGenDB(db)
}
\arguments{
  \item{db}{Handle of the database to be closed} 
}
\details{
	This function closes an existing open database whose handle is provided by the input argument.
}
\seealso{
\code{\link{createGenDB}},
\code{\link{openGenDB}}
}
\references{

}
\examples{
SQLite()
drv<-dbDriver("SQLite");
db <- openGenDB("16S.db")
closeGenDB(db)
}
\keyword{manip}
