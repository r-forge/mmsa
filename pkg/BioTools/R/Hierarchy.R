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
GenClass16S <- function(domain=NA, phylum=NA, class=NA, order=NA, 
	family=NA, genus=NA, species=NA, strain=NA, id=NA) {

    c(Domain=domain, Phylum=phylum, Class=class,
	    Order=order, Family=family, Genus=genus, Species=species,
	    Strain=strain, id=id)
}

### get classification info
getTaxonomyNames <- function(db) dbListFields(db$db, "classification")

getRank <- function(db, rank="Phylum", whereRank=NULL, whereName=NULL, 
	table="sequences", all=FALSE, partialMatch = TRUE, 
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
	    partialMatch=FALSE
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
	    .getWhere(db, whereRank, whereName, partialMatch), " GROUP BY classification.", cols, 
                     " ORDER BY count(classification.",cols, ") desc", sep="" )
    else
	statement <- paste("SELECT ", distinct, " classification.",cols,
	    getIDs,
	    " FROM ",table," t INNER JOIN classification ON t.id=classification.id ", 
	    .getWhere(db, whereRank, whereName, partialMatch), sep="")

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
	partialMatch = TRUE, removeUnknown=FALSE) 
getRank(db, rank="id", whereRank, whereName, table, all=TRUE, partialMatch, removeUnknown)


getHierarchy <- function(db, rank, name, drop=TRUE, partialMatch=TRUE){
    hierarchy <- getTaxonomyNames(db)

    if(missing(name)|| missing(rank)) stop("no name and/or rank given!")

    .getHierarchy <- function(db, rank, name) {
	#convert rank to number, eg: kingdom=1, phylum=2
	rankNum <- which(tolower(hierarchy)==tolower(.pmatchRank(db,rank)))

	cl <- sapply(1:rankNum, FUN=function(i) 
		as.character(getRank(db, rank=hierarchy[i], 
				whereRank=rank, whereName=name,
				partialMatch=FALSE)))

	m <- matrix(NA, nrow=1, ncol=length(hierarchy))
	m[1:rankNum] <- unlist(cl)[1:rankNum]
	m
    }

    
    ### find all matching names
    name <- unlist(lapply(name, FUN=function(n) 
		    	getRank(db, rank=rank, whereRank=rank, whereName=n, 
			    partialMatch=partialMatch)))
    
    if(length(name) >0) { ### handle multiple names
    m <- t(sapply(name, FUN=function(x) 
		    .getHierarchy(db, rank, x)))
    colnames(m) <- getTaxonomyNames(db)
    rownames(m) <- NULL
    }else{
        coln <- getTaxonomyNames(db)
	m <- data.frame(rep(list(character(0)), length(coln)))
	colnames(m) <- coln
    }

    if(drop) m <- drop(m)
    m
}


### helper
.pmatchRank <- function(col, rank, numeric=FALSE) {
    fields <- dbListFields(col$db, "classification")
    if(is.null(rank)) m <- 1
    else {
	m <- pmatch(tolower(rank), tolower(fields))
	if(any(is.na(m))) stop("Rank not found!")
    }

    if(numeric) m else fields[m]
}

.getWhere <- function(db, rank, name, partialMatch=TRUE, removeUnknownSpecies=FALSE) {
    if(partialMatch) exact <- "%" else exact <- ""

    if(is.null(rank) && is.null(name)) where <- ""
    else if (length(name) <=1)  where <- paste("WHERE classification.'", .pmatchRank(db, rank), 
		"' LIKE '", name, exact, "'", sep='')
	
	#more than one names are provided
	else if (length(name) > 1) 
 	{
		rankExact <- .pmatchRank(db, rank)


		### find all matching names
		name <- unlist(lapply(name, FUN=function(n) 
				getRank(db, rank=rank, whereRank=rank, 
					whereName=n, 
					partialMatch=partialMatch)))
    
		where <- paste("WHERE classification.'", rankExact, 
			"' IN ('", paste(name,collapse="','"), "')", sep='')
	
	}
	if (removeUnknownSpecies)
		where <- paste(where, " AND classification.species NOT LIKE 'Unknown%' ")
	where
}

### helper if we want to use a formula interface for where and name
.formulaTowhere <- function(formula=NULL) {
  if(is.null(formula)) {
    rank<-NULL
    name<- NULL
  }else{
  f <- as.character(eval(formula))
  if(f[[1L]]!="~" || length(f)!=3) stop("Incorrect formula!")
 
    rank <- f[[2L]]
    name <- unlist(strsplit(f[[3L]], split="\\s*\\+\\s*"))
  }
  
  list(rank=rank, name=name)
}
