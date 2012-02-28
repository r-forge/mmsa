### implement simrank similarity defined in
### Santis et al, Simrank: Rapid and sensitive general-purpose k-mer 
### search tool, BMC Ecology 2011, 11:11

### x and y are vectors of sequences
simRank <- function(x, y, k=3, alphabet = "acgt") {
    x.kmers <- sapply(x, FUN=function(s) 
	    count(s2c(s), word=k, alphabet=s2c(alphabet))) >0
    y.kmers <- sapply(y, FUN=function(s) 
	    count(s2c(s), word=k, alphabet=s2c(alphabet))) >0

    shared <- crossprod(x.kmers, y.kmers)
    min <- outer(colSums(x.kmers), colSums(y.kmers), function(x, y) { 
		sapply(1:length(x), FUN = function(i) min(x[i], y[i]))
	    })

    shared/min
}

