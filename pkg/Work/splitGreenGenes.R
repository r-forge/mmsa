library("MMSA")
library("seqinr")

.parseAnnotationFull <- function(annot)
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

f<- "current_GREENGENES_gg16S_unaligned.fasta"
cat("doing", f, "\n")
    # read the fasta formatted file 
    if(!is(try(sequences <- read.fasta(f)), "try-error")){
      cat("Length of file is ",length(sequences),"\n")
      for(i in 1:length(sequences))
      {
        annot <- getAnnot(sequences[[i]])
        annot <- sub(">", "", annot)
		cl<-.parseAnnotationFull(annot)
		#phylum
		seqClass<-cl["p__"]
		filename<-paste("phylum/",seqClass,".fasta",sep="")				
        if(!file.exists(filename))
          file.create(filename)
        write(getAnnot(sequences[[i]]),file=filename,append=T)
        write(toupper(getSequence(sequences[[i]],as.string=T)),file=filename,append=T)
		#class
		seqClass<-cl["c__"]
		filename<-paste("class/",seqClass,".fasta",sep="")				
        if(!file.exists(filename))
          file.create(filename)
        write(getAnnot(sequences[[i]]),file=filename,append=T)
        write(toupper(getSequence(sequences[[i]],as.string=T)),file=filename,append=T)
      	#order
		seqOrder<-cl["o__"]
		filename<-paste("order/",seqOrder,".fasta",sep="")				
        if(!file.exists(filename))
          file.create(filename)
        write(getAnnot(sequences[[i]]),file=filename,append=T)
        write(toupper(getSequence(sequences[[i]],as.string=T)),file=filename,append=T)
		#family
		seqFamily<-cl["f__"]
		filename<-paste("family/",seqFamily,".fasta",sep="")				
        if(!file.exists(filename))
          file.create(filename)
        write(getAnnot(sequences[[i]]),file=filename,append=T)
        write(toupper(getSequence(sequences[[i]],as.string=T)),file=filename,append=T)
		#genus
		seqFamily<-cl["g__"]
		filename<-paste("genus/",seqFamily,".fasta",sep="")				
        if(!file.exists(filename))
          file.create(filename)
        write(getAnnot(sequences[[i]]),file=filename,append=T)
        write(toupper(getSequence(sequences[[i]],as.string=T)),file=filename,append=T)
		#species  
		seqSpecies<-cl["s__"]
		seqSpecies<-gsub("/","",seqSpecies)
		filename<-paste("species/",seqSpecies,".fasta",sep="")				
        if(!file.exists(filename))
          file.create(filename)
        write(getAnnot(sequences[[i]]),file=filename,append=T)
        write(toupper(getSequence(sequences[[i]],as.string=T)),file=filename,append=T)
	
	  }
    }





