### implement simrank similarity defined in
### Santis et al, Simrank: Rapid and sensitive general-purpose k-mer 
### search tool, BMC Ecology 2011, 11:11

### x and y are vectors of sequences
simRank <- function(x, k=7) {
    x <- DNAStringSet(x)
    
    x.kmer <- oligonucleotideFrequency(x, k) >0

    shared <- tcrossprod(x.kmer)
    sum.kmer <- rowSums(x.kmer)
    
    min <- outer(sum.kmer, sum.kmer, function(x, y) pmin(x,y)) 

    as.simil(shared/min)
}

