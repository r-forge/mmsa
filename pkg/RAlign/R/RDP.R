#This function requires either Qiime or MacQiime to be installed.

classifyRDP <- function(x, param=NULL){

    ## get temp files and change working directory
    wd <- tempdir()
    dir <- getwd()
    temp_file <- basename(tempfile(tmpdir = wd))
    on.exit({
		file.remove(Sys.glob(paste(temp_file, "*", sep="")))
		setwd(dir)
	    })
    setwd(wd)

    infile <- paste(temp_file, ".fasta", sep="")
    outfile <- paste(temp_file, "_tax_assignments.txt", sep="")

    write.XStringSet(x, infile, append=FALSE, format="fasta")

    system(paste(.findExecuable(c("qiime", "macqiime")), 
		    "assign_taxonomy", "-i", infile, "-o ." , "-m rdp"))	

    cl_tab <- read.table(outfile, sep="\t") 
    
    ### FIXME: is 10 good?
    cl <- as.data.frame(t(sapply(strsplit(as.character(cl_tab[,2]), ";"), 
		    function(x) { length(x) <- 10; x })))
    
    rownames(cl) <- cl_tab[,1]
    colnames(cl) <- paste("level", 1:10, sep="_")
    cl <- cbind(cl, prop=cl_tab[,3])
    
    ### rdp does not seem to produce the output in the same order!
    cl <- cl[names(x),]
    cl
}


### FIXME: this needs to go to QuasiAlign
### Maybe the database stuff of quasialign should go to Ralign?
#checkRDPOutput <- function(db,RDPfile="rdp/classify_tax_assignments.txt")
#{
#	if (!file.exists(RDPfile))
#		stop("File not found")
#	correctFile <- "rdp/correct.txt"
#	incorrectFile <- "rdp/incorrect.txt"
#	if (file.exists(correctFile))
#		file.remove(correctFile)
#	if (file.exists(incorrectFile))
#		file.remove(incorrectFile)
#	cat("Id","RDP Classification","Correct Classification","\n",sep="\t",append=TRUE,file=incorrectFile)
#	cat("Id","RDP Classification","Correct Classification","\n",sep="\t",append=TRUE,file=correctFile)
#	correct <- 0
#	rdp<- read.table(RDPfile, stringsAsFactors=FALSE)
#	ids <- rdp[,1]
#	h<- getHierarchy(db,rank="id",name=ids,partialMatch=FALSE)
#	for(i in 1:nrow(rdp))
#	{
#		correctClassification <- gsub(" \\(class\\)","",paste(h[i,1:6],collapse=";"))
#		if (length(grep(correctClassification,rdp[i,2])) >0)
#			{
#				correct <- correct + 1
#				cat(rdp[i,1],rdp[i,2],correctClassification,"\n",sep="\t",append=TRUE,file=correctFile)
#			}
#		else
#				cat(rdp[i,1],rdp[i,2],correctClassification,"\n",sep="\t",append=TRUE,file=incorrectFile)
#	}
#
#	cat("Total sequences = ",nrow(rdp),"Correct classification = ",correct," Percentage  =",(correct/nrow(rdp))*100,"\n")
#}
