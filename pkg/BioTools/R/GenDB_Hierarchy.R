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


GenClass16S <- function(Domain=NA, Phylum=NA, Class=NA, Order=NA, 
                        Family=NA, Genus=NA, Species=NA, Strain=NA, Id=NA) {
  
  ### prevent recycling
  params <- as.list(environment())
  l <- max(sapply(params, length))
  params <- lapply(params, "length<-", l)
  return(as.data.frame(do.call(cbind, params)))
}
  

### this is a helper to create sequences with only the id as names
Annotation_Id <- function(annotation, decode) {
  if(decode) { ### decode metadata
    stop("Not a decoder")
  }else{ ### recreate meta data   
    ret <- annotation[,"Id"]
  }
  
  ret
}




### get classification info
getTaxonomyNames <- function(db) dbListFields(db$db, "classification")

getRank <- function(db, rank="Phylum", whereRank=NULL, whereName=NULL, 
                    table="sequences", all=FALSE, 
                    count=FALSE, removeUnknown=FALSE) {
  
  
  fields <- getTaxonomyNames(db)
  cols <- paste("[", fields[.pmatchRank(db, rank, 
                                        numeric=TRUE)],"]", sep='')
  rankPosition <- .pmatchRank(db, rank, numeric=TRUE) 
  if (!is.null(whereRank))
  {
    whereRankPosition <- .pmatchRank(db, whereRank, numeric=TRUE) 
    
    ### for id we always return all
    if (whereRankPosition == .pmatchRank(db,"id",numeric=TRUE))
    {
      all=TRUE
    }
  }
  else whereRankPosition = 0
  
  if(all) {
    distinct <- ""
    getIDs <- paste(", classification.", .pmatchRank(db, whereRank) ,sep="")
  }else{
    distinct <- "DISTINCT"
    getIDs <- ""
  }
  
  
  if(count)
    statement <- paste("SELECT ", distinct, " classification.",cols,
                       " , count(classification.", cols, ") AS count FROM ", table," t INNER JOIN classification ON t.id=classification.id  ", 
                       .getWhere(db, whereRank, whereName), " GROUP BY classification.", cols, 
                       " ORDER BY count(classification.",cols, ") desc", sep="" )
  else
    statement <- paste("SELECT ", distinct, " classification.",cols,
                       getIDs,
                       " FROM ",table," t INNER JOIN classification ON t.id=classification.id ", 
                       .getWhere(db, whereRank, whereName), sep="")
  
  ret <- dbGetQuery(db$db, statement = statement)
  
  
  if(removeUnknown)
    ret <- ret[ret[,1]!="unknown" & ret[,1]!="NA",, drop=FALSE]
  
  
  if(nrow(ret)<1) warning("No matching name in rank found!")
  
  #count
  if(count) return(structure(ret[,2], names=ret[,1]))
  
  #unique values
  if(!all) return(ret[,1])    
  
  #ids
  if(rank=="id") return(ret[,1])
  
  #match order
  if(nrow(ret)==length(whereName)) {
    # if the number of results in the same as the number of names (reorder) 
    ret <- ret[match(whereName,ret[,2]),,drop=FALSE]
    return(structure(as.factor(ret[,1]), names=ret[,2]))
  }
  
  #default
  return(as.factor(ret[,1]))
  
}

getIDs <- function(db, whereRank=NULL, whereName=NULL, table="sequences",  
                   removeUnknown=FALSE) { 
  getRank(db, rank="id", whereRank, whereName, table, all=TRUE, removeUnknown)
}

### FIXME: recode for speed
getHierarchy <- function(db, rank, name){
  hierarchy <- getTaxonomyNames(db)
  
  ### shortcut if rank is a DNAStringSet produced with get_sequences
  if(is(rank, "DNAStringSet")) {
    name <- names(rank)
    rank <- "id"
  }
  
  if(missing(name)|| missing(rank)) stop("no name and/or rank given!")
  
  statement <- paste("SELECT * FROM classification", 
                     .getWhere(db, rank, name)) 
  res <- dbGetQuery(db$db, statement = statement)
  
  if(nrow(res)<1) warning("No matching hierarchy found!")

  ### make sure the oder is correct if rank is id
  if (!is.null(rank) && rank=="id" && nrow(res)==length(name)) {
    res<-res[match(name,res$Id),]  
  }
  
  
  res
}


### helper
.pmatchRank <- function(db, rank, numeric=FALSE) {
  fields <- getTaxonomyNames(db)
  
  if(is.null(rank)) m <- 1
  else {
    m <- pmatch(tolower(rank), tolower(fields))
    if(any(is.na(m))) stop("Rank not found!")
  }
  
  if(numeric) m else fields[m]
}

.getWhere <- function(db, rank, name, removeUnknownSpecies=FALSE) {
  
  if(is.null(rank) || is.null(name)) {
    if(removeUnknownSpecies) return("classification.species NOT LIKE 'Unknown%'")
    else return("")
  }
  
  if(length(rank)!=1) stop("Exactly one rank has to be given!")
  rank <- .pmatchRank(db, rank)
  
  if(length(name)<=1) {
    ### partial match using %?
    if(length(grep("%", name))) exact <- " LIKE " else exact <- "="
    where <- paste("WHERE classification.'", rank, "'", exact, "'", name, "'", 
                   sep='')
  }else{
    ### partial match using %?
    if(length(grep("%", name))) {
      name <- unlist(lapply(name, FUN=function(n) 
        getRank(db, rank=rank, whereRank=rank, whereName=n)))
    }
    where <- paste("WHERE classification.'", rank, "' IN ('", 
                   paste(name, collapse="','"), "')", sep='')
  }
  
  if (removeUnknownSpecies)
    where <- paste(where, "AND classification.species NOT LIKE 'Unknown%'")
  
  where
}

