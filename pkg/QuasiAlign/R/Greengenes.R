
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


GreengenesMetaDataReader <- function(annotation) {

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

    c(cl, org_name,id)
}


