NSV_counter <- function(x, window=50, overlap=0, word=3, 
        offset = 0, last_window=FALSE) {

    x <- getSequence(x)
    
    l <- as.integer(length(x)/(window-overlap)) -1L
    start <- (window-overlap)*(0:l) + 1
    end <- start + window - 1

    ## fixme: test offset
    if(offset != 0) {
        start <- start + offset
        end <- end + offset

        start <- start[start <= length(x)]
        end <- end[end <= length(x)]
    }

    if(last_window) {
        if(tail(end,1) < length(x)) {
            start <- c(start, tail(end,1)+1)
            end <- c(end, length(x))
        }
    }

    t(sapply(1:length(start), FUN=function(i) 
                    count(x[start[i]:end[i]], word=word)))
}


NSV_list_counter <- function(x, window=50, overlap=0, word=3, 
        offset=0, last_window=FALSE) lapply(x, NSV_counter, 
        window=window, overlap=overlap, word=word, offset=offset,
        last_window=last_window)


NSV_make_stream <- function(cnt, use_ss=TRUE, ss_val = NA) {
    ## start state
    ss <- NULL
    if(use_ss) ss <- rep(ss_val, ncol(cnt[[1]]))

    stream <- matrix(NA, ncol= ncol(cnt[[1]]), nrow=0)
    colnames(stream) <- colnames(cnt[[1]])

    for(i in 1:length(cnt)) stream <- rbind(stream, ss, cnt[[i]])

    stream
}

