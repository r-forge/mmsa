library(Biostrings)

## install clustalw (is a Ubuntu package) 

## example
filepath <- "test2_org.FASTA"
s <- read.DNAStringSet(filepath, format="fasta",
	nrec=-1L, skip=0L, use.names=TRUE)
al <- clustal(s)
al


clustal <- function(x) {
    ## get temp files and change working directory
    wd <- tempdir()
    dir <- getwd()
    
    on.exit(setwd(dir))
    setwd(wd)

    temp_file <- basename(tempfile(tmpdir = wd))

    write.XStringSet(s, temp_file, append=FALSE, format="fasta")

    ## call clustalw (needs to be installed and in the path!)
    system(paste("clustalw", temp_file))

    read.DNAMultipleAlignment(paste(temp_file, ".aln", sep=""), 
	    format="clustal")
}
