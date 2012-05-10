\name{Greengenes}
\alias{GenClass16S}
\alias{GenClass16S_Greengenes}
\alias{GreengenesMetaDataReader}

\title{Manage GenClass16S objects}
\description{
	These functions create objects of class GenClass16S which contain classification hierarchy.
}
\usage{
	GenClass16S_Greengenes(kingdom=NA, phylum=NA, class=NA, order=NA, 
	    family=NA, genus=NA, species=NA, otu=NA, org_name=NA, id=NA) 
	GreengenesMetaDataReader(annotation)	
}
\arguments{
	\item{kingdom}{Name of the kingdom to which the organism belongs. } 
	\item{phylum}{Name of the phylum to which the organism belongs. } 
	\item{class}{Name of the class to which the organism belongs. } 
	\item{order}{Name of the order to which the organism belongs. } 
	\item{family}{Name of the family to which the organism belongs. } 
	\item{genus}{Name of the genus to which the organism belongs. } 
	\item{species}{Name of the species to which the organism belongs. } 
	\item{otu}{Name of the otu to which the organism belongs. } 
	\item{org_name}{Name of the organism. } 
	\item{id}{ID of the sequence from Greengenes. } 
	\item{annotation}{description line for a sequence in FASTA format.}
}
\details{
\code{GenClass16S_Greengenes} creates a classification object for a 
	sequence. The creator is also used for defining
	the classification structure when a new GenDB is created.

\code{GreengenesMetaDataReader} converts the description line found
in Greengenes FASTA files into a \code{GenClass16S_Greengenes} object 
suitable to be stored in a GenDB.
The reader is used in \code{addSequences}.
}
\seealso{
\code{\link{addSequences}},
\code{\link{createGenDB}}
}
\examples{
### example of how to create a classification object
GenClass16S_Greengenes(
	kingdom="Bacteria", phylum ="Firmicutes", class="Clostridia", 
	order="Clostridiales", family="Peptococcaceae",
	genus="Desulfotomaculum", species="Unclassified", otu="2143", 
	org_name="Desulfotomaculum thermocisternum str. ST90", id="13651")
}
