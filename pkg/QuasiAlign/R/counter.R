#start
.counter <- function(x, window=100, overlap=0, word=3, 
	last_window=FALSE, allOffsets=FALSE) {
    
    x <- DNAString(x)
    if (!allOffsets)
	{
	#returns the sequence as a vector  

    ret <- matrix(NA, ncol=4^word, nrow=0)
    
    if(length(x) < window) {
		warning("Sequence is shorter than window size!")
		return(matrix(nrow=0, ncol=0))
   		}
	l <- as.integer(length(x)/(window-overlap)) -1L
	#start and end positions as vectors
	start <- (window-overlap)*(0:l) + 1
	end <- start + window - 1

	if(last_window) {
	    if(tail(end,1) < length(x))	{
		start <- c(start, tail(end,1)+1)
		end <- c(end, length(x))
	    }
	}

	mat <- t(sapply(1:length(start), FUN=function(i) 
			oligonucleotideFrequency(DNAString(x,
					start=start[i],
					nchar=end[i]-start[i]+1), 
				word)))
	ret <- rbind(ret,mat)
	}	
	else if (allOffsets)
	{
		freq<-oligonucleotideFrequency(DNAString(x,start=1,nchar=window),word)
		fmatrix <- matrix(nrow=(length(x)-window+1),ncol=4^word)
		fmatrix[1,]<-as.matrix(t(freq))
		colnames(fmatrix) <- colnames(as.matrix(t(freq)))
		start<-2
		end <- window+1
		while(end <=length(x))
		{
				next_mer<-DNAString(x,start=end-(word-1),nchar=word)
				first_mer<-DNAString(x,start=start-1,nchar=word)
				#update
				freq[as.character(first_mer)]<-freq[as.character(first_mer)]-1
				freq[as.character(next_mer)]<-freq[as.character(next_mer)]+1
				#f<-rbind(f,t(freq)[,1:(4^word)])
				fmatrix[start,]<-as.matrix(t(freq)[,1:(4^word)])
				start<-start+1
				end<-end+1
		}
		rearrangedIdx<-unlist(lapply(1:window,FUN=function(j){c(seq(j, length(x)-window+1, by=100),NA)}))
		ret <- fmatrix[rearrangedIdx,]	
	}
	ret

}

#   if (allOffsets)
#	end <- window
#    else
#	end <- 1
#
##    for(seqOffset in 1:end)
#    {   
#	y <- DNAString(x, start=seqOffset)    
#	if (length(y) < window) break;
#	l <- as.integer(length(y)/(window-overlap)) -1L
#
#	#start and end positions as vectors
#	start <- (window-overlap)*(0:l) + 1
#	end <- start + window - 1
#
#	if(last_window) {
#	    if(tail(end,1) < length(x))	{
#		start <- c(start, tail(end,1)+1)
#		end <- c(end, length(x))
#	    }
#	}
#
#	mat <- t(sapply(1:length(start), FUN=function(i) 
#			oligonucleotideFrequency(DNAString(y,
#					start=start[i],
#					nchar=end[i]-start[i]+1), 
#				word)))
#	
#	if (seqOffset>1) ret <- rbind(ret, matrix(NA,nrow=1,ncol=4^word))
#	ret <- rbind(ret, mat)  
#    }
#    ret
#}

#end


# count individual sequences
.count_sequences <- function(x, window=100, overlap=0, 
	word=3, last_window=FALSE) 
lapply(x, .counter, window=window, overlap=overlap, word=word, 
	last_window=last_window)


## take several sequences and append them with the start state
## symbol
.make_stream <- function(cnt, use_ss=TRUE, ss_val = NA) {
    ## start state
    ss <- NULL
    if(use_ss) ss <- rep(ss_val, ncol(cnt[[1]]))

    stream <- matrix(NA, ncol= ncol(cnt[[1]]), nrow=0)
    colnames(stream) <- colnames(cnt[[1]])

    for(i in 1:length(cnt)) stream <- base::rbind(stream, ss, cnt[[i]])

    stream
}

