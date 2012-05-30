

clustal <- function(x) {
    ## get temp files and change working directory
    wd <- tempdir()
    dir <- getwd()
    temp_file <- basename(tempfile(tmpdir = wd))
    
    on.exit({file.remove(temp_file); setwd(dir)})
    setwd(wd)


    write.XStringSet(x, temp_file, append=FALSE, format="fasta")

    ## call clustalw (needs to be installed and in the path!)
    system(paste(Sys.which("clustalw"), temp_file))

    read.DNAMultipleAlignment(paste(temp_file, ".aln", sep=""), 
	    format="clustal")
}
