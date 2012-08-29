## classification hierarchy for 16S
GenClass16S <- function(domain=NA, phylum=NA, class=NA, order=NA, 
	family=NA, genus=NA, species=NA, strain=NA) {

    c(Domain=domain, Phylum=phylum, Class=class,
	    Order=order, Family=family, Genus=genus, Species=species,
	    Strain=strain)
}

### get classification info
getTaxonomyNames <- function(db) {
    cl <- dbListFields(db$db, "classification")
    cl <- head(cl, length(cl))   ### remove data
    cl
}

getRank <- function(db, rank=NULL, whereRank=NULL, whereName=NULL, 
	all=FALSE, partialMatch = TRUE, count=FALSE, removeUnknown=FALSE) {

    fields <- getTaxonomyNames(db)
    cols <- paste("[", fields[.pmatchRank(db, rank, 
		    numeric=TRUE)],"]", sep='')
	
	if(all) distinct <- "" else  distinct <- "DISTINCT"

	if(count)
		statement <- paste("SELECT ", distinct, " classification.",cols,
		    " , count(", cols, ") AS count FROM classification ", 
		    .getWhere(db, whereRank, whereName, partialMatch), " GROUP BY ", cols, " ORDER BY count(",cols, ") desc" )
	else
		statement <- paste("SELECT ", distinct, " classification.",cols,
		    " FROM classification ", 
		    .getWhere(db, whereRank, whereName, partialMatch))
	
    ret <- dbGetQuery(db$db, 
	    statement = statement)
	if(removeUnknown)
		if (length(which(ret[1,]=="unknown")) > 0)
			ret <- ret[-which(ret[1,]=="unknown"),]
	ret
}

getHierarchy <- function(db, rank, name, drop=TRUE, partialMatch=TRUE){
    hierarchy <- getTaxonomyNames(db)
    
    .getHierarchy <- function(db, rank, name) {
		#convert rank to number, eg: kingdom=1, phylum=2
		rankNum <- which(tolower(hierarchy)==tolower(.pmatchRank(db,rank)))
		
		cl <- sapply(1:rankNum, FUN=function(i) 
			getRank(db, rank=hierarchy[i], whereRank=rank, whereName=name,
				partialMatch=FALSE)[,1])

		m <- matrix(NA, nrow=1, ncol=length(hierarchy))
		m[1:rankNum] <- unlist(cl)[1:rankNum]
		m
    }

    ### find all matching names
    name <- unlist(lapply(name, FUN=function(n) 
		    	getRank(db, rank=rank, whereRank=rank, whereName=n, 
			    partialMatch=partialMatch)[,1]))
    

    if(length(name) < 1) stop("No match found!")

    ### handle multiple names
    m <- t(sapply(name, FUN=function(x) 
		    .getHierarchy(db, rank, x)))
    colnames(m) <- getTaxonomyNames(db)
    rownames(m) <- NULL

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

.getWhere <- function(col, rank, name, partialMatch=TRUE) {
    if(partialMatch) exact <- "%" else exact <- ""

    if(is.null(rank) && is.null(name)) where <- ""
    else if (length(name) <=1)  where <- paste("WHERE classification.'", .pmatchRank(col, rank), 
		"' LIKE '", name, exact, "'", sep='')
	#more than one names are provided
	else if (length(name) > 1)  where <- paste("WHERE classification.'", .pmatchRank(col, rank), 
		"' IN (", paste(name,collapse=","), ")", sep='')
	
	where
}


