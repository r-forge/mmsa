modelStatesPlot <- function (model, states, ylab=TRUE, ...)
{
		window <- as.integer(model$window)
		#xrange is the max number of windows from all the sequences
		xrange <- max(sapply(model$clusterInfo, FUN=function(x) length(x)))*window
		#yrange is the number of sequences
		yrange <- model$nSequences
		plot(NA, xlim=c(0,xrange),ylim=c(1,max(1.05*yrange,2)), axes=FALSE, xlab="nucleotide positions", ylab="sequences", bty="l", las=1, ...)
		#bty="l"
		#labels for x axis
		axis(1, at=seq(0,xrange,by=window), labels=seq(0,xrange,by=window), las=2)
		if(ylab)
			axis(2, at=seq(0,yrange,by=5),labels =seq(0,yrange,by=5),las=1)
		else
			axis(2,las=1)
		abline(h=c(1:yrange), col="grey", lty="dotted")
		colors <- c("red")
		vertColors <- c("blue")
		colorIdx <- 0	
	for(modelState in states)
	{
		#get meta info from clustering details of the model
		cd <- getModelDetails(model, modelState)
		#get the sequences which are part of the state 
		sequences <- cd[,1]
		#get the indexes of the sequences
		#sequenceInd <- which(names(model$clusterInfo) %in% sequences)
		sequenceInd <- vector()
		for(i in 1:length(sequences)) { sequenceInd[i] <- which(names(model$clusterInfo)==sequences[i]) }
		#get the segments
		segments <- cd[,2]
		#colorIdx <- colorIdx + 1
		#color <- colors[(colorIdx+1)%%2 + 1 ]
		#vertColor <- vertColors[(colorIdx+1)%%2 + 1]
		#vertical lines	
		#lines(as.numeric(segments),sequenceInd, col=vertColors, lwd=2)
		#lines(as.numeric(segments)+1,sequenceInd, col=vertColors, lwd=2)
		#lines(as.numeric(segments)+0.5,sequenceInd, lty=2)
		lines(((segments-1)*window+1)+window/2,sequenceInd, lty=2)
		#horizontal lines
		for(i in 1:length(sequenceInd))
				lines(.segmentToSequenceNumbers(model, segments[i]),c(sequenceInd[i],sequenceInd[i]), col=colors, lwd=2)
		#put the text for the state number at the top
		text(((segments[length(segments)]-1)*window+1)+window/2, as.numeric(max(sequenceInd)),modelState, pos=3, adj=c(0,0), xpd=TRUE, col="black")
	}
	hyper<-list(c(69,99), c(137,242), c(433,497), c(576,682), c(822,879), c(986,1043), c(1117,1173), c(1243,1294),c(1435,1465))
	for(i in 1:length(hyper))
	{	
		#lines(hyper[[i]],c(max(sequenceInd)+1,max(sequenceInd)+1),col="blue",lwd=2)
		lines(hyper[[i]],c(max(1.05*yrange,2),max(1.05*yrange,2)),col="blue",lwd=2)
		hyperRegion <- paste("V",i,sep="")
		#text(mean(hyper[[i]]),max(sequenceInd)+1, hyperRegion ,pos=3, adj=c(0,0), xpd=TRUE, col="black")
		text(mean(hyper[[i]]),max(1.05*yrange,2), hyperRegion ,pos=3, adj=c(0,0), xpd=TRUE, col="black")
			
	}

}

.segmentToSequenceNumbers<- function(model, segment)
{
	window <- as.integer(model$window)	
	start <- (segment-1)*window + 1
	end <- start + window -1 
	#snum <- list()
	#for(i in 1:length(start)) {
	#	snum[[i]] <- c(start[i],end[i])
	#}
	snum <- c(start,end)
	names(snum) <- c("start","end")
	return(snum)	
}

compareSequences <- function(model, sequences, ...)
{
		#xrange is the max number of windows from all the sequences
		xrange <- max(sapply(model$clusterInfo, FUN=function(x) length(x)))
		#yrange is the number of sequences
		yrange <- model$nSequences
		plot(NA, xlim=c(1,xrange), ylim=c(1,yrange), axes=FALSE, xlab="segments", ylab="sequences", ...)
		#labels for x axis
		axis(1, at=1:xrange, labels=c(1:xrange))
		#labels for y axis
		axis(2, at=1:yrange, labels=names(model$clusterInfo))
		#horizontal lines
		abline(h=c(1:yrange), col="grey", lty="dotted")

		if (length(sequences) < 2) stop("Need to compare at least 2 sequences")
	
		ci <- model$clusterInfo
		common <- intersect( ci[[sequences[1]]], ci[[sequences[2]]] )
		
		if(length(sequences) >= 3)
			for(i in seq(3,length(sequences))) common <- intersect(common, ci[[sequences[i]]])	
		
		if (length(common)==0) print("No common states found between those sequences")
		
		colors <- c("red")
		vertColors <- c("blue")
		colorIdx <- 0	
		for(modelState in common)
		{
			#get meta info from clustering details of the model
			cd <- getModelDetails(model, modelState)
			#get the sequences 
			modelSequences <- cd[,1]
			#get the indexes of the sequences in the model
			sequenceInd <- which(names(model$clusterInfo) %in% modelSequences)
			#find missing indexes
			missing <- which(!(sequenceInd %in% sequences))	
			#filter down to just the required
			sequenceInd <- intersect(sequenceInd, sequences)
			#get the segments
			segments <- as.numeric(cd[,2])
			if (length(missing) > 0 )
				segments <- segments[-missing]
			colorIdx <- colorIdx + 1
			color <- colors[(colorIdx+1)%%2 + 1 ]
			vertColor <- vertColors[(colorIdx+1)%%2 + 1]
			#vertical blue lines
			#lines(segments,sequenceInd, col=vertColors, lwd=2)
			#lines(segments+1,sequenceInd, col=vertColors, lwd=2)
			lines(segments+0.5,sequenceInd, lty=2)
			#horizontal lines
			for(i in 1:length(sequenceInd))
					lines(c(segments[i],segments[i]+1),c(sequenceInd[i],sequenceInd[i]), col=colors, lwd=2)
			#put the text for the state number at the top
			text(as.numeric(segments[length(segments)])+0.5 ,as.numeric(max(sequenceInd)), modelState, pos=3, xpd=TRUE, col="black")
			#if( (max(sequenceInd) != 1) )	
			#	text(as.numeric(segments[length(segments)])+0.35,as.numeric(max(sequenceInd))-0.10,modelState, col="black")
			#else if (max(sequenceInd) == 1)
			#	text(as.numeric(segments[length(segments)])+0.35,as.numeric(max(sequenceInd))+0.10,modelState, col="black")
		}


}

findLargestCommon <- function(model, limit=NULL)
{
	largestCommon <- list()
	ci <- model$clusterInfo
	if(!is.null(limit)) 
		limit <-min(limit,model$nSequences)
	else
		limit <- model$nSequences

	for(r in seq(2,limit))
	{
		#take r combinations from the names vector 
		comb <- combn(names(ci), r)
		maxCommon <- 0
		maxList <-vector()
		for(i in 1:ncol(comb))
		{
			#take 2 at a time
			common <- intersect(unlist(ci[comb[1,i]]),unlist(ci[comb[2,i]]))
			if (nrow(comb) >= 3)
			{
				for(j in 3:nrow(comb))
				{	
					common <- intersect(common, unlist(ci[comb[j,i]]))
				}
			}
			
			if (length(common) > maxCommon) 
				{
					maxCommon <- length(common)
					maxList <- as.vector(comb[,i])
				}

		}
			#cat("size = ",r," has max common = ",maxCommon,"\n")
			#print(maxList)			
			largestCommon[[r]]<-which(names(ci) %in% maxList)
			
	}
	return(largestCommon)
}

