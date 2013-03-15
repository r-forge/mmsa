
createRDPTraining <- function(db, rank, name, ids)
{
	tax <- getHierarchy(db,rank="id",name=ids)
	#tax is a matrix
	outTrainTax <- vector()
	if (!file.exists("rdp/taxonomy"))
		dir.create("rdp/taxonomy", recursive=TRUE)
	if (!file.exists("rdp/sequences"))
		dir.create("rdp/sequences", recursive=TRUE)
	for(i in 1:nrow(tax))
	{
		header =""
		#seq = ""
		id = tax[,"Id"][i]
		seq = getSequences(db,rank="id",name=id,partialMatch=FALSE)
		header<- paste(id,"\t",sep="")
		header <- paste(header,"Root",sep="")
		for(j in 2:6)
		{
			header <- paste(header, ";" ,substr(tolower(names(tax[i,j])),1,1),"__",tax[i,j],sep="")
		}
		header <- gsub("\\(class\\)","",header)
		header <- gsub(" ","",header)		
		
		outTrainTax[i] <- header
		write.XStringSet(seq,filepath=paste("rdp/sequences/",rank,name,".fasta",sep=""), append=TRUE)
	}
	write(outTrainTax, file=paste("rdp/taxonomy/",rank,name,".txt",sep=""), append=TRUE)
}

