
classifyRDP <- function(x, confidence=.8, property=NULL, java_args="-Xmx1g"){

    ## check 
    if(Sys.getenv("RDP_JAR_PATH") =="") stop("Environment variable 'RDP_JAR_PATH needs to be set!'")

    ## get temp files and change working directory
    #wd <- tempdir()
    #dir <- getwd()
    #temp_file <- basename(tempfile(tmpdir = wd))
    #on.exit({
	#	file.remove(Sys.glob(paste(temp_file, "*", sep="")))
	#	setwd(dir)
	#    })
    #setwd(wd)

    #infile <- paste(temp_file, ".fasta", sep="")
    #outfile <- paste(temp_file, "_tax_assignments.txt", sep="")
   
	infile <-"query.fasta" 
    outfile <-"outFile.out"
	 ## property?
    if(!is.null(property)) property <- paste("-t", property)
    else property <- ""

    write.XStringSet(x, infile, append=FALSE)

    system(paste("java", java_args, "-jar", Sys.getenv("RDP_JAR_PATH"), 
		    property, "-q", infile, "-o", outfile),
	    ignore.stdout=TRUE, ignore.stderr=TRUE)

    ## read and parse rdp output
    cl_tab <- read.table(outfile, sep="\t") 
    
    ## remove empty columns
    cl_tab <- cl_tab[!sapply(cl_tab, FUN=function(x) all(is.na(x)))]
    
    seq_names <- cl_tab[,1] ## sequence names are in first column
    
    i <- seq(2, ncol(cl_tab), by=3) ## 3 columns for each tax. level
    
    ## get classification
    cl <- cl_tab[,i]	
    dimnames(cl) <- list(seq_names, as.matrix(cl_tab[1,i+1])[1,])

    ## get confidence
    conf <- as.matrix(cl_tab[,i+2])
    dimnames(conf) <- list(seq_names, as.matrix(cl_tab[1,i+1])[1,])
    
    if(confidence>0) cl[conf < confidence] <- NA
    
    attr(cl, "confidence") <- conf    
    cl
}

trainRDP <- function(sequences, java_args="-Xmx1g", outDir=".") 
{
    if(Sys.getenv("RDP_JAR_PATH") =="") stop("Environment variable 'RDP_JAR_PATH needs to be set!'")
	
	write.XStringSet(sequences,"train.fasta")
	l<-strsplit(names(sequences),"Root;")
	annot<-sapply(l,FUN=function(x) x[2])
	h<-matrix(ncol=6,nrow=0)
	colnames(h) <-c("Kingdom","Phylum","Class","Order","Family","Genus")
	for(i in 1:length(annot)) {h<-rbind(h,unlist(strsplit(annot[i],";"))[1:6])}
	m<-matrix(ncol=5,nrow=0)
	#first row of the file
	f<-"0*Root*-1*0*rootrank"
	m<-rbind(m,unlist(strsplit(f,split="\\*")))
	taxNames <- colnames(h)
	for(i in 1:nrow(h))
	{
		for(j in 1:ncol(h))	
		{
			taxId <- nrow(m)
			taxonName <- h[i,j]
			if (j==1)
				parentTaxId <-0
			else
				parentTaxId <- previousTaxId
			depth <- j
			rank <- colnames(h)[j]
		
			#search if already there
			if (length(which(m[,2]==taxonName & m[,5]==rank)) == 0)
			{
				str <- paste(taxId,taxonName,parentTaxId,depth,rank,sep="*")
				m<-rbind(m,unlist(strsplit(str,split="\\*")))
				previousTaxId <- taxId
			}	
			else if (length(which(m[,2]==taxonName & m[,5]==rank)) > 0)
			{
				row <- which(m[,2]==taxonName & m[,5]==rank)
				#cat("Row =",row,"\n")
				#print(m[row,1])
				previousTaxId <- m[row,1]
			}
			#end seach
		}
	}
	out<-apply(m,MARGIN=1,FUN=function(x) paste(x,collapse="*"))
	write(out, file="train.txt")
	#create parsed training files
	system(paste("java", java_args, "-cp", Sys.getenv("RDP_JAR_PATH"),"edu/msu/cme/rdp/classifier/train/ClassifierTraineeMaker train.txt train.fasta 1 version1 test ", outDir),
	    ignore.stdout=TRUE, ignore.stderr=TRUE)
	file.copy(system.file("examples/rRNAClassifier.properties",package="RAlign"),outDir)
	return(file.path(outDir,"rRNAClassifier.properties"))


}

