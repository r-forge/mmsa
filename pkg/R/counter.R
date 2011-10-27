library("seqinr")
#x 
#word => size of mer
counter <- function(x, window=50, overlap=0, word=3, last_window=FALSE) {
  #returns the sequence as a vector  
  x <- getSequence(x)    
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
#count function is part of seqinr package - counts occurence of "word" sized mers
    t(sapply(1:length(start), FUN=function(i) 
            count(x[start[i]:end[i]], word=word)))
}
# count individual sequences
count_sequences <- function(x, window=50, overlap=0, word=2, last_window=FALSE) 
lapply(x, counter, window=window, overlap=overlap, word=word, 
    last_window=last_window)


make_stream <- function(cnt, use_ss=TRUE, ss_val = NA) {
    ## start state
    ss <- NULL
	if(use_ss) ss <- rep(ss_val, ncol(cnt[[1]]))

    stream <- matrix(NA, ncol= ncol(cnt[[1]]), nrow=0)
    colnames(stream) <- colnames(cnt[[1]])

	for(i in 1:length(cnt)) stream <- rbind(stream, ss, cnt[[i]])
    
    stream
}

create_stream <- function(file, window=100, overlap=0, 
    word=3, last_window=FALSE, max=0) {
    ## read file
    sequences <- read.fasta(file) 

    if(max!=0 && length(sequences)>max) sequences <- sequences[1:max]
 dnafile <- system.file("sequences/someORF.fsa", package = "seqinr")
    ## count
    cnt <- count_sequences(sequences, 
        window=window, overlap=overlap, word=word, last_window=last_window)

    ## report number of windows
    cat(file, "\n")
    cat("number of sequences", length(sequences) ,"\n")
    cat("window length distr.", "\n")
    print(summary(sapply(cnt, nrow)))

    ## put together and add start states
    stream <- make_stream(cnt, use_ss = TRUE, ss_val=NA)
    attr(stream, "n") <- length(sequences)

    stream
}

