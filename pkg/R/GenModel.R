
genModel <- function(db, rank=NULL, name=NULL, table, limit=-1, 
	measure="Kullback", threshold=0.10, plus_one=TRUE) {
	
	### FIXME: check in metadata if it is NSV

	emm <- EMM(measure=measure,threshold=threshold)	
	d<-getSequences(db, rank, name, table, limit=limit)
	
	for(i in 1:length(d$NSV))
	{
		#sequence<-decodeSequence(d$NSV[i])
		sequence<- d$NSV[i]		
		if(plus_one) sequence <- sequence +data.frame(1)
		build(emm,sequence)
		reset(emm)
	}	

	# add name attribute
	rank <- .pmatchRank(db, rank)
	name <- d[[rank]][1]
	op <- paste(rank,": ", name, " - ", nrow(d), " sequences" , sep = '')
	attr(emm, "name") <- op	
	emm    
	
}

plot.GenModel<-function(emm,...)
{
	plot(emm, main=attr(emm, "name"),...)
	
}

