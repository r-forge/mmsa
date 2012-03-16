
## classification hierarchy for 16S
GenClass16S_Greengenes <- function(kingdom=NA, phylum=NA, class=NA, order=NA, 
	family=NA, genus=NA, species=NA, otu=NA,org_name=NA,id=NA) {

    c(Kingdom=kingdom, Phylum=phylum, Class=class,
	    Order=order, Family=family, Genus=genus, Species=species,
	    Otu=otu,Org_name=org_name,Id=id)
}





addSequencesGreengenes <- function(db, file, verbose=FALSE) {

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

	c(cl, org_name,id)
    }
	#end helper function
    ok <- 0
    fail <- 0
    total <- 0

    dbBeginTransaction(db$db)
    #start
	f <- read.DNAStringSet(file)
	for(i in 1:length(f)) {
		annot<- names(f)[i]
		cl <- .parseAnnotation(annot)
		org_name<-cl[length(cl)]
		cl <- paste("'",cl,"'", sep='', collapse=', ') 
		#print(cl)
		dat<- f[[i]]
		dat<- tolower(as.character(dat[1:length(dat)]))
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

    #close(f)
    
    cat("Read", ok+fail, "entries. Added", ok , "entries.\n")

}



