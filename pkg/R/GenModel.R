
genModel <- function(db, rank=NULL, name=NULL, table, limit=-1, 
	measure="Kullback", threshold=0.10, plus_one=TRUE) {
	
	### FIXME: check in metadata if it is NSV

	emm <- EMM(measure=measure,threshold=threshold)	
	d<-getSequences(db, rank, name, table, limit=limit)
	
	for(i in 1:length(d))
	{
		#sequence<-decodeSequence(d$NSV[i])
		sequence<- d[i]		
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

plot.genModel<-function(x, ...)
{
	plot(x, main=attr(emm, "name"),...)
	
}

