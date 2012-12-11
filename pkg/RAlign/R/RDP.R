#######################################################################
# RAlign - Interfaces to several sequence alignment and classification tools
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

classifyRDP <- function(x, confidence=.8, classifier=NULL, java_args="-Xmx1g"){

    ## check 
    if(Sys.getenv("RDP_JAR_PATH") =="") stop("Environment variable 'RDP_JAR_PATH needs to be set!'")

   # get temp files and change working directory
   wd <- tempdir()
   dir <- getwd()
   temp_file <- tempfile(tmpdir = wd)
   #temp_file <- "train"
   on.exit({
				file.remove(Sys.glob(paste(temp_file, "*", sep="")))
				setwd(dir)
			})

	#setwd(wd)
    infile <- paste(temp_file, ".fasta", sep="")
    outfile <- paste(temp_file, "_tax_assignments.txt", sep="")
	 ## property?
    if(!is.null(classifier)) property <- paste("-t", file.path(classifier,"rRNAClassifier.properties"))
    else property <- ""

    writeXStringSet(x, infile, append=FALSE)
	if (system(paste("java", java_args, "-jar", Sys.getenv("RDP_JAR_PATH"), 
		    property, "-q", infile, "-o", outfile),
	    ignore.stdout=TRUE, ignore.stderr=TRUE)) 
		stop("Error executing jar: ",Sys.getenv("RDP_JAR_PATH"))

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

trainRDP <- function(x, java_args="-Xmx1g", classifier="classifier") 
{
    if(Sys.getenv("RDP_JAR_PATH") =="") stop("Environment variable 'RDP_JAR_PATH needs to be set!'")
	if (length(list.files(classifier)) > 0) stop("Classifier directory should be empty")
	if (!file.exists("classifier/")) dir.create("classifier")
	
	writeXStringSet(x,file.path(classifier,"train.fasta"))
	l<-strsplit(names(x),"Root;")
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
	write(out, file=file.path(classifier,"train.txt"))
	#create parsed training files
	system(paste("java", java_args, "-cp", Sys.getenv("RDP_JAR_PATH"),"edu/msu/cme/rdp/classifier/train/ClassifierTraineeMaker ",file.path(classifier,"train.txt"), file.path(classifier,"train.fasta")," 1 version1 test ", classifier),
	    ignore.stdout=TRUE, ignore.stderr=TRUE)
	file.copy(system.file("examples/rRNAClassifier.properties",package="RAlign"),classifier)
	
	invisible(classifier)


}

findAccuracy <- function(actual, predicted, rank)
{
	rank <- colnames(actual)[grep(rank, colnames(actual), ignore.case=TRUE)]
	actual <- factor(actual[,rank])
    predicted <- factor(predicted[,rank])

	commonLevels <- sort(unique(c(levels(actual), levels(predicted))))
	actual <- factor(actual, levels = commonLevels)
	predicted <- factor(predicted, levels = commonLevels)

	table(actual,predicted, dnn=list("actual","predicted"))
}
