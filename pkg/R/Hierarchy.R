## classification hierarchy for 16S
GenClass16S <- function(domain=NA, phylum=NA, class=NA, order=NA, 
	family=NA, genus=NA, species=NA, strain=NA) {

    c(Domain=domain, Phylum=phylum, Class=class,
	    Order=order, Family=family, Genus=genus, Species=species,
	    Strain=strain)
}

### get classification info
getClassification <- function(db) {
    cl <- dbListFields(db$db, "classification")
    cl <- head(cl, length(cl))   ### remove data
    cl
}

getRank <- function(db, rank=NULL, whereRank=NULL, whereName=NULL, 
	all=FALSE, partialMatch = TRUE) {
    fields <- getClassification(db)
    cols <- paste("[", fields[.pmatchRank(db, rank, 
		    numeric=TRUE)],"]", sep='')

    if(all) distinct <- "" else  distinct <- "DISTINCT"

    dbGetQuery(db$db, 
	    statement = paste("SELECT", distinct, "classification.",cols,
		    " FROM classification ", 
		    .getWhere(db, whereRank, whereName, partialMatch), " ORDER BY ",cols))
}

getHierarchy <- function(db, rank, name, drop=TRUE, partialMatch=TRUE){
    hierarchy <- getClassification(db)
    
    .getHierarchy <- function(db, rank, name) {
		#convert rank to number, eg: kingdom=1, phylum=2
		rankNum <- which(tolower(hierarchy)==tolower(.pmatchRank(db,rank)))
		
		cl <- sapply(1:rankNum, FUN=function(i) 
			getRank(db, rank=hierarchy[i], whereRank=rank, whereName=name,
				partialMatch=FALSE))

		m <- matrix(NA, nrow=1, ncol=length(hierarchy))
		m[1:rankNum] <- unlist(cl)[1:rankNum]
		m
    }

    ### find all matching names
    name <- unlist(lapply(name, FUN=function(n) 
		    getRank(db, rank=rank, whereRank=rank, whereName=n, 
			    partialMatch=partialMatch)))
    

    if(length(name) < 1) stop("No match found!")

    ### handle multiple names
    m <- t(sapply(name, FUN=function(x) 
		    .getHierarchy(db, rank, x)))
    colnames(m) <- getClassification(db)
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
    else where <- paste("WHERE classification.'", .pmatchRank(col, rank), 
		"' LIKE '", name, exact, "'", sep='')
    where
}


