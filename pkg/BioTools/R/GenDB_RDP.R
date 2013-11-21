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
GenClass16S_RDP <- GenClass16S_Greengenes

Annotation_RDP <- function(annotation, decode=TRUE) {
  
  if(decode) { ### decode metadata
    
    ann <- strsplit(annotation, ";")
    
    ret <- matrix(NA_character_, ncol=length(GenClass16S_RDP()), nrow=length(ann))
    for(i in 1:length(ann)) ret[i,1:length(ann[[i]])] <- ann[[i]] 
    
    ret[,1] <- sub(" Root$", "", ret[,1])
    
    ret <- GenClass16S_RDP(
      Kingdom=as.character(ret[,2]), 
      Phylum=as.character(ret[,3]), 
      Class=as.character(ret[,4]), 
      Order=as.character(ret[,5]),                  
      Family=as.character(ret[,6]), 
      Genus=as.character(ret[,7]), 
      Species=NA_character_, 
      Otu=NA_character_,
      Org_name=NA_character_,
      Id=as.character(ret[,1])
    )
    
  }else{ ### recreate meta data   
    h <- annotation[,1:6]
    h <- apply(h,MARGIN=1,FUN=function(x) paste(x,collapse=";"))
    h <- gsub(";unknown","", h)
    h <- gsub(" \\(class\\)","", h)
    
    ret <- paste(annotation[, "Id"], " ", "Root;", h, sep="")
    
  }
  
  ret
}


