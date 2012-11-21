#This function requires either Qiime or MacQiime to be installed.

classifyRDP <- function(sequences=NULL, filepath=NULL)
{
	#time <-gsub("/","",Sys.time())
	#time <-gsub(" ","",time)
	filename <- paste("rdp/classify.fasta",sep='')
	if (!file.exists("rdp/"))
    	dir.create("rdp/", recursive=TRUE)
	if(!is.null(sequences))
    	write.XStringSet(sequences, filepath=filename)	
	else if(!is.null(filepath))
		file.copy(filepath,filename)
	if (Sys.which("macqiime") !="")
		exec <- "macqiime"
	else if (Sys.which("qiime") !="")
		exec <- "qiime"
	else
	{
		print("unable to find macqiime/qiime for RDP")
		exec <- NULL
	}
	if(!is.null(exec))
	{
		command <- paste("assign_taxonomy.py -i ",filename," -o rdp/ -m rdp")
		system(exec, input=command)
	}
}

checkRDPOutput <- function(db,RDPfile="rdp/classify_tax_assignments.txt")
{
	if (!file.exists(RDPfile))
		stop("File not found")
	correctFile <- "rdp/correct.txt"
	incorrectFile <- "rdp/incorrect.txt"
	if (file.exists(correctFile))
		file.remove(correctFile)
	if (file.exists(incorrectFile))
		file.remove(incorrectFile)
	cat("Id","RDP Classification","Correct Classification","\n",sep="\t",append=TRUE,file=incorrectFile)
	cat("Id","RDP Classification","Correct Classification","\n",sep="\t",append=TRUE,file=correctFile)
	correct <- 0
	rdp<- read.table(RDPfile, stringsAsFactors=FALSE)
	ids <- rdp[,1]
	h<- getHierarchy(db,rank="id",name=ids,partialMatch=FALSE)
	for(i in 1:nrow(rdp))
	{
		correctClassification <- gsub(" \\(class\\)","",paste(h[i,1:6],collapse=";"))
		if (length(grep(correctClassification,rdp[i,2])) >0)
			{
				correct <- correct + 1
				cat(rdp[i,1],rdp[i,2],correctClassification,"\n",sep="\t",append=TRUE,file=correctFile)
			}
		else
				cat(rdp[i,1],rdp[i,2],correctClassification,"\n",sep="\t",append=TRUE,file=incorrectFile)
	}

	cat("Total sequences = ",nrow(rdp),"Correct classification = ",correct," Percentage  =",(correct/nrow(rdp))*100,"\n")
}
