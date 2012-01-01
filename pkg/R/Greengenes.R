
## classification hierarchy for 16S
GenClass16S_Greengenes <- function(kingdom=NA, phylum=NA, class=NA, order=NA, 
	family=NA, genus=NA, species=NA, otu=NA,org_name=NA) {

    c(kingdom=kingdom, phylum=phylum, class=class,
	    order=order, family=family, genus=genus, species=species,
	    otu=otu,org_name=org_name)
}





addSequencesGreengenes <- function(db, file,verbose=FALSE) {

    #helper function
    .parseAnnotation <- function(annot)
    {
	
	fields <- c("k__", "p__", "c__", "o__", "f__", "g__", "s__", "otu_") 
	## add out and org_name
	
	# remove leading ">"
	annot <- sub(">", "", annot)
	# split at "k__" 
	tmp <- strsplit(annot, " *k__")[[1]]
	org_name <- gsub("'","",tmp[1])	
	tmp <- strsplit(paste('k__',tmp[2], sep=''), '; *')[[1]]
	cl <- sapply(fields, FUN=function(f) {
		    val <- grep(f, tmp, value=TRUE)
		    val <- sub('^.__', '', val)
		    if(length(val) ==0) val <- "unknown"
		    val
		})

	c(cl, org_name)
    }
	#end helper function
    ok <- 0
    fail <- 0
    total <- 0

    dbBeginTransaction(db$db)
    f <- file(file)
    open(f)

    while(TRUE)  {
	
	annot <- ''
	while(substr(annot,1,1) != '>') {	
	    annot <- readLines(f,1)
	    if(length(annot) ==0) break
	}
	dat <- tolower(readLines(f,1))
	
	if(length(annot) ==0) break
	cl <- .parseAnnotation(annot)
	org_name<-cl[length(cl)]
	
	cl <- paste("'",cl,"'", sep='', collapse=', ') 
	


	tr <- try(dbSendQuery(db$db,          
			statement = paste("INSERT INTO classification VALUES(", 
				cl,  ")", sep='')), silent=FALSE)

	
	
	tr <- try(dbSendQuery(db$db,          
			statement = paste("INSERT INTO sequences VALUES('", 
				org_name, "','", dat,  "')", sep='')), 
		silent=FALSE)


	if(!is(tr, "try-error")) ok <- ok+1
	else fail <- fail+1
	total <- total+1
	if(verbose)
		{	
			if(total%%100 == 0) cat("Read", total, "entries (ok:", ok, 
				"/ fail:", fail,")\n")
		}
    }
    dbCommit(db$db)

    close(f)
    
    cat("Read", ok+fail, "entries. Added", ok , "entries.\n")

}



