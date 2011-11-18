
## classification hierarchy for 16S
GenClass16S_Greengenes <- function(kingdom=NA, phylum=NA, class=NA, order=NA, 
	family=NA, genus=NA, species=NA, otu=NA,org_name=NA) {

    c(kingdom=kingdom, phylum=phylum, class=class,
	    order=order, family=family, genus=genus, species=species,
	    otu=otu,org_name=org_name)
}


addSequencesGreengenes <- function(collection, file, verb=FALSE) {

    #helper function
    .readSequence <- function(currSequence)
    {
	#start processing
	annot <- getAnnot(currSequence[[1]])
	# remove leading ">"
	annot <- sub(">", "", annot)
	#split based on ;
	x<-gregexpr(";",annot)
	#start1 for splitting the headers
	classification=list()
	classification[[1]]=c("kingdom","phylum","class","order","family","genus","species","otu")
	classification[[2]]=c("k__","p__","c__","o__","f__","g__","s__","otu_")
	classification[[3]]=""
	for(j in 1:length(classification[[1]]))
	    #for(j in 1:8)
	{
	    location <-regexpr(text=annot,pattern=classification[[2]][j])
	    if(location[1]!=-1)
	    {
		y<-which(x[[1]]>location[1])
		if (length(y) == 0)
		    y<-nchar(annot)
		end<- x[[1]][y[1]]
		if (is.na(end))
		    end <- nchar(annot)
		if (j<8)
		    classification[[3]][j] <- substr(annot,location[1]+3,end-1)
		else if (j==8)
		    classification[[3]][j] <- substr(annot,location[1],end-1)
		#remove _
		sub("_","",classification[[3]][j])
		if (is.na(classification[[3]][j]) || classification[[3]][j]=="")
		{
		    classification[[3]][j]<- 'UNKNOWN'            
		}            
	    }
	    else
		classification[[3]][j]='UNKNOWN'
	    cat(paste(classification[[3]][j],"\t"),file="summary.txt",append=T)


	} #for (j in 1:8)
	cat("\n",file="summary.txt",append=T)
	org_name<-strsplit(annot,split=" k__")[[1]][1]
	org_name<-sub(";","",org_name)
	org_name <- sub(" ","",org_name)
	kingdom = classification[[3]][1]
	phylum = classification[[3]][2]
	class = classification[[3]][3]
	order = classification[[3]][4]
	family = classification[[3]][5]
	genus = classification[[3]][6]
	species = classification[[3]][7]
	otu = classification[[3]][8]
	#end1
	#make a GenClass16S_Greengenes class object
	gen16class<- GenClass16S_Greengenes(kingdom,phylum,class,order,family,genus,species,otu,org_name)
	sequence<- getSequence(currSequence[[1]],as.string=TRUE)[[1]]
	return(list(sequence=sequence,classification=gen16class))


	#end processing

    }


    ok <- 0
    fail <- 0

    if(file.info(file)[1,"isdir"]) files <- dir(file, full.names=TRUE)
    else files <- file

    dbBeginTransaction(collection$db)
    for(f in files){
	if(verb) cat("Processing:", f, "\n")
	if(!is(try(sequences <- read.fasta(f)), "try-error")){
	    for(i in 1:length(sequences)){

		sequence <-.readSequence(sequences[i])

		cl <- paste("'",sequence$classification,"'", sep='', 
			collapse=', ')
		org_name <-paste("'",sequence$classification[length(sequence$classification)],"'",sep='')

		## FIXME: NSVs?
		dat <- sequence$sequence

		## Insert into DB
		#tr <- try(dbSendQuery(collection$db,          
		#		statement = paste("INSERT INTO ",
		#			collection$collection, " VALUES(", 
		#			cl, ", '", dat, "')", sep='')), silent=FALSE)
		
		tr <- try(dbSendQuery(collection$db,          
				statement = paste("INSERT INTO ",
					collection$collection, " VALUES(", 
					cl,  ")", sep='')), silent=FALSE)
		
		tr <- try(dbSendQuery(collection$db,          
				statement = paste("INSERT INTO ",
					paste(collection$collection,"Seq",sep=""), " VALUES(", 
					org_name, ",'", dat,  "')", sep='')), silent=FALSE)

		if(!is(tr, "try-error")) ok <- ok+1
		else fail <- fail+1


		}
	}
    }
    dbCommit(collection$db)

    cat("Read", ok+fail, "entries. Added", ok , "entries to", 
	    collection$collection,".\n")

}




addSequencesGreengenes_large <- function(collection, file) {

    #helper function
    .parseAnnotation <- function(annot)
    {
	
	fields <- c("k__", "p__", "c__", "o__", "f__", "g__", "s__", "otu_") 
	## add out and org_name
	
	# remove leading ">"
	annot <- sub(">", "", annot)
	# split at "k__" 
	tmp <- strsplit(annot, " *k__")[[1]]
	org_name <- tmp[1]
	tmp <- strsplit(paste('k__',tmp[2], sep=''), '; *')[[1]]
	cl <- sapply(fields, FUN=function(f) {
		    val <- grep(f, tmp, value=TRUE)
		    val <- sub('^.__', '', val)
		    if(length(val) ==0) val <- "unknown"
		    val
		})

	c(cl, org_name)
    }

    ok <- 0
    fail <- 0
    total <- 0

    dbBeginTransaction(collection$db)
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
	cl <- paste("'",cl,"'", sep='', collapse=', ') 


	## Insert into DB
	tr <- try(dbSendQuery(collection$db,          
			statement = paste("INSERT INTO ",
				collection$collection, " VALUES(", 
				cl, ", '", dat, "')", sep='')), silent=TRUE)

	if(!is(tr, "try-error")) ok <- ok+1
	else fail <- fail+1
	total <- total+1

	if(total%%100 == 0) cat("Read", total, "entries (ok:", ok, 
		"/ fail:", fail,")\n")

    }
    dbCommit(collection$db)

    close(f)
    
    cat("Read", ok+fail, "entries. Added", ok , "entries to", 
	    collection$collection,"\n")

}



