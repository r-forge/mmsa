\name{GenDB_Greengenes}
\alias{GenDB_Greengenes}
\alias{GenClass16S}
\alias{GenClass16S_Greengenes}


\title{Manage GenClass16S objects}
\description{
	These functions create objects of class GenClass16S which contain classification hierarchy.
}
\usage{
	GenClass16S_Greengenes(Kingdom=NA, Phylum=NA, Class=NA, Order=NA, 
	    Family=NA, Genus=NA, Species=NA, Otu=NA, Org_name=NA, Id=NA) 

}
\arguments{
	\item{Kingdom}{Name of the kingdom to which the organism belongs. } 
	\item{Phylum}{Name of the phylum to which the organism belongs. } 
	\item{Class}{Name of the class to which the organism belongs. } 
	\item{Order}{Name of the order to which the organism belongs. } 
	\item{Family}{Name of the family to which the organism belongs. } 
	\item{Genus}{Name of the genus to which the organism belongs. } 
	\item{Species}{Name of the species to which the organism belongs. } 
	\item{Otu}{Name of the otu to which the organism belongs. } 
	\item{Org_name}{Name of the organism. } 
	\item{Id}{ID of the sequence from Greengenes. } 
}
\details{
\code{GenClass16S_Greengenes} creates a classification object for a 
	sequence. The creator is also used for defining
	the classification structure when a new GenDB is created.

\code{Annotation_Greengenes} is used to decode and create Greengenes-format
annotations.
}
\seealso{
  \code{\link{GenDB}},
  \code{\link{GenDB_Hierarchy}},
  \code{\link{Sequences}}
}
\examples{
### example of how to create a classification object
GenClass16S_Greengenes(
	Kingdom="Bacteria", Phylum ="Firmicutes", Class="Clostridia", 
	Order="Clostridiales", Family="Peptococcaceae",
	Genus="Desulfotomaculum", Species="Unclassified", Otu="2143", 
	Org_name="Desulfotomaculum thermocisternum str. ST90", Id="13651")
}
