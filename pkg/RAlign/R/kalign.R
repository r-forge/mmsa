kalign <- function(x, param=NULL) {
    
    ## get temp files and change working directory
    wd <- tempdir()
    dir <- getwd()
    temp_file <- basename(tempfile(tmpdir = wd))
    on.exit({
		file.remove(Sys.glob(paste(temp_file, ".*", sep=""))) 
		setwd(dir)
	    })
    setwd(wd)

    infile <- paste(temp_file, ".in", sep="")
    outfile <- paste(temp_file, ".aln", sep="")
    
    write.XStringSet(x, infile, append=FALSE, format="fasta")

    ## call clustalw (needs to be installed and in the path!)
    system(paste(.findExecuable("kalign"), infile, outfile, "-f fasta", param))

    read.DNAMultipleAlignment(outfile, format="fasta")
}

kalign_help <- function() {
    system(paste(.findExecuable("kalign"), "-h"))
}

