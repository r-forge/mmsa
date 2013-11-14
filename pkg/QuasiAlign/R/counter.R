#start
.counter <- function(x, window=100L, overlap=0L, word=3L, 
                     last_window=FALSE, allOffsets=FALSE) {
  
  if(!is(x, "DNAString")) stop("x needs to be a DNAString!")
  
  if (!allOffsets){
    if(length(x) < window) {
      warning("Sequence is shorter than window size!")
      return(matrix(NA_integer_, nrow=0L, ncol=0L))
    }
    
    l <- as.integer(length(x)/(window-overlap)) -1L
    #start and end positions as vectors
    start <- (window-overlap)*(0:l) + 1L
    end <- start + window - 1L
    
    if(last_window) {
      if(tail(end,1L) < length(x))	{
        start <- c(start, tail(end,1L)+1L)
        end <- c(end, length(x))
      }
    }
    
    mat <- t(sapply(1:length(start), FUN=function(i) 
      oligonucleotideFrequency(DNAString(x,
                                         start=start[i],
                                         nchar=end[i]-start[i]+1L), 
                               word)))
    return(mat)
  }

  ### allOffset is TRUE or 1
  if((is.logical(allOffsets) && allOffsets) || allOffsets==1) { 
    
    if(last_window) warning("last_window is ignored for allOffsets")
    
    ### create matrix and get count for initial segment
    freq <- oligonucleotideFrequency(DNAString(x, start=1L, nchar=window), word)
    
    ### subsetting arrays of characters is much faster than for DNAString
    xc <- BioTools:::s2c(as.character(x)) ### sequence as a character array
    xc[!(xc %in% c("G","A","T", "C"))] <- NA ### ignor odd characters
    n <- length(xc)
    
    ### create count matrix
    mat <- matrix(NA_integer_, nrow=n+1L, ncol=4^word) #n-window+1 + window NAs
    colnames (mat) <- names(freq)
    mat[1,] <- freq
    
    start <- 2L
    end <- window+1L
    
    ### the counts are stored in the order 1st pass NA, 2nd pass, NA, etc.
    ### but the counting is done all 1st segments, all 2nd segments, etc.
    ### number of segments for each run plus NA separator row
    pos <- c(1L, cumsum(floor(n:(n-window+2L)/window)+1L)+1L)
    
    
    rearrangedIdx <- unlist(lapply(1:window, FUN=function(j) 
      seq.int(j, n-window+1L, by=window)))
    
    while(end <= n){
      remove_word <- BioTools:::c2s(xc[(start-1L):(start+word-2L)])
      add_word <- BioTools:::c2s(xc[(end-word+1L):end])
      
      if(nchar(remove_word)==word) freq[remove_word] <- freq[remove_word]-1L
      if(nchar(add_word)==word) freq[add_word] <- freq[add_word]+1L
      
      mat[pos[start%%window]+floor(start/window),] <- freq
      
      start<-start+1L
      end<-end+1L
    }
  
    return(mat)
  }
  
  ### allOffset is a number
  if(last_window) warning("last_window is ignored for allOffsets")
  
  mat <- sapply(seq(1L, window, by=allOffsets), FUN=function(i) {
      s <- DNAString(x, start=i)
      rbind(.counter(s, window=window, overlap=overlap, word=word, 
               last_window=FALSE, allOffsets=FALSE), NA_integer_)
      })
  
  mat <- do.call(rbind, mat)
  return(mat)
  
}


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

