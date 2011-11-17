
genModel <- function(collection,rank=NULL,name=NULL, n=-1, measure="Kullback", 
	threshold=0.10, plus_one=TRUE) {
	#if (collection$collection ! ="nsv")
	#	stop("Not in NSV format")
	emm <- EMM(measure=measure,threshold=threshold)	
	d<-getSequences(nsv, rank, name)
	for(i in 1:length(d$sequence))
	{
		sequence<-decodeSequence(d$sequence[i])
		if(plus_one) sequence <- sequence +1
		build(emm,sequence)
		reset(emm)
	}	

	# add name attribute
	rank <- .pmatchRank(collection, rank)
	name <- d[[rank]][1]
	op <- paste(rank,": ", name, " - ", nrow(d), " sequences" , sep = '')
	attr(emm, "name") <- op	
	emm    
	
}

plot.GenModel<-function(emm,...)
{
	plot(emm, main=attr(emm, "name"),...)
	
}

