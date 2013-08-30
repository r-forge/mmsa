#######################################################################
# BioTools - Interfaces to several sequence alignment and 
# classification tools
# Copyright (C) 2012 Michael Hahsler and Anurag Nagar
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

## classification hierarchy for 16S
GenClass16S_Greengenes <- function(kingdom=NA, phylum=NA, class=NA, order=NA, 
	family=NA, genus=NA, species=NA, otu=NA,org_name=NA,id=NA) {

    c(	    Kingdom=as.character(kingdom), 
	    Phylum=as.character(phylum), 
	    Class=as.character(class),
	    Order=as.character(order), 
	    Family=as.character(family), 
	    Genus=as.character(genus), 
	    Species=as.character(species),
	    Otu=as.character(otu),
	    Org_name=as.character(org_name),
	    Id=as.character(id)
	    )
}


Annotation_Greengenes <- function(annotation, decode=TRUE) {

  if(decode) { ### decode metadata
    fields <- c("k__", "p__", "c__", "o__", "f__", "g__", "s__", "otu_") 
    ## add out and org_name

    # remove leading ">"
    annotation <- sub(">", "", annotation)
    # split at "k__" 
    tmp <- strsplit(annotation, " *k__")[[1]]
    org_name <- gsub("'","",tmp[1])
    id<-strsplit(org_name," ")[[1]][1]
    #org_name<- trimSpace(sub(id,"",org_name))
    org_name<- gsub(" ","",(sub(id,"",org_name)))
    tmp <- strsplit(paste('k__',tmp[2], sep=''), '; *')[[1]]
    cl <- sapply(fields, FUN=function(f) {
		val <- grep(f, tmp, value=TRUE)
		val <- sub('^.__', '', val)
		if(length(val) ==0) val <- "unknown"
		val
	    })

    ret <- c(cl, org_name,id)
  
  }else{ ### recreate meta data   
    ret <- paste(
      ">", annotation[,"Id"], " ", annotation[,"Org_name"],
      " k__", annotation[,"Kingdom"], 
      " p__", annotation[,"Phylum"],            
      " c__", annotation[,"Class"], 
      " o__", annotation[,"Order"], 
      " f__", annotation[,"Family"],
      " g__", annotation[,"Genus"], 
      " s__", annotation[,"Species"], 
      " otu_", annotation[,"Otu"], 
      sep="")
  }
  
  ret
  }


