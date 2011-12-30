\name{GenClass16S}
\alias{GenClass16S}
\alias{GenClass16S_Greengenes}
\docType{data}
\title{
	Manage GenClass16S objects
}
\description{
	These functions create objects of class GenClass16S which contain classification hierarchy.
}
\usage{
	GenClass16S_Greengenes(kingdom, phylum, class, order, family, genus, species, otu, org_name) 
}
\arguments{
	\item{kingdom}{Name of the kingdom to which the organism belongs. OPTIONAL} 
	\item{phylum}{Name of the phylum to which the organism belongs. OPTIONAL} 
	\item{class}{Name of the class to which the organism belongs. OPTIONAL} 
	\item{order}{Name of the order to which the organism belongs. OPTIONAL} 
	\item{family}{Name of the family to which the organism belongs. OPTIONAL} 
	\item{genus}{Name of the genus to which the organism belongs. OPTIONAL} 
	\item{species}{Name of the species to which the organism belongs. OPTIONAL} 
	\item{otu}{Name of the otu to which the organism belongs. OPTIONAL} 
	\item{org_name}{Name of the organism. OPTIONAL} 
}
\details{
	The function GenClass16S_Greengenes creates a classification object for the gene sequence. 
}
\seealso{
\code{\link{addSequencesGreengenes}},
\code{\link{nSequences}},
\code{\link{getSequences}},
\code{\link{getClassification}}
}
\examples{

	GenClass16S_Greengenes(kingdom="Bacteria", phylum ="Firmicutes", class="Clostridia", order="Clostridiales", family="Peptococcaceae"
	,genus="Desulfotomaculum", species="Unclassified", otu="2143", org_name="Desulfotomaculum thermocisternum str. ST90" )

}
\keyword{datasets}
