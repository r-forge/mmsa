library(QuasiAlign)

context("Counter")

seqs <- readDNAStringSet( system.file("examples/phylums/Proteobacteria.fasta", 
		                              package="QuasiAlign"))

s <- seqs[1:2]
NSV <- createNSVSet(s)

### check positions 1 and 10
expect_identical(NSV[[1]][1,], oligonucleotideFrequency(
  DNAString(s[[1]], start=1, nchar=100), width=3))

expect_identical(NSV[[1]][10,], oligonucleotideFrequency(
  DNAString(s[[1]], start=100*(10-1)+1, nchar=100), width=3))

expect_identical(NSV[[2]][10,], oligonucleotideFrequency(
  DNAString(s[[2]], start=100*(10-1)+1, nchar=100), width=3))

### all Offsets
NSV2 <- createNSVSet(s, allOffsets=TRUE)
n_segments <- nrow(NSV[[1]]) 

### first n segments should be identical
expect_identical(NSV2[[1]][1:n_segments,], NSV[[1]][1:n_segments,])
### check for NA separator
expect_equivalent(NSV2[[1]][n_segments+1,], rep(NA_integer_, 64))

### check new pass start
expect_identical(NSV2[[1]][n_segments+2,], oligonucleotideFrequency(
  DNAString(s[[1]], start=2, nchar=100), width=3))

### should end wirh NA separator
expect_equivalent(NSV2[[1]][nrow(NSV2[[1]]),], rep(NA_integer_, 64))


### check second sequence
expect_identical(NSV2[[2]][1,], NSV[[2]][1,])



