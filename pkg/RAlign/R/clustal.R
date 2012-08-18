clustal <- function(x, param=NULL) {
    
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
    reader <- if(is(x, "RNAStringSet")) read.RNAMultipleAlignment
	else if(is(x, "DNAStringSet")) read.DNAMultipleAlignment
	else if(is(x, "AAStringSet")) read.AAMultipleAlignment
	else stop("Unknown sequence type!")


    write.XStringSet(x, infile, append=FALSE, format="fasta")

    ## call clustalw (needs to be installed and in the path!)
    system(paste(.findExecuable("clustalw"), infile, param))

    reader(outfile, format="clustal")
}

clustal_help <- function() {
    system(paste(.findExecuable("clustalw"), "-help"))
}


clustal_profile <- function(x, y, param=NULL) {
    
    if(is(y, "DNAMultipleAlignment")) {
	profileprofile <- TRUE
	param2 <- "-profile"
    } else {
	profileprofile <- FALSE
	param2 <- "-sequence"
    }

    ## get temp files and change working directory
    wd <- tempdir()
    dir <- getwd()
    temp_file <- basename(tempfile(tmpdir = wd))
    on.exit({
		file.remove(Sys.glob(paste(temp_file, ".*", sep=""))) 
		setwd(dir)
	    })
    setwd(wd)

    prof1 <- paste(temp_file, ".in1", sep="")
    prof2 <- paste(temp_file, ".in2", sep="")
    outfile <- paste(temp_file, ".aln", sep="")
    
    write.phylip(x, prof1)
    if(profileprofile) write.phylip(y, prof2)
    else write.XStringSet(y, prof2, append=FALSE, format="fasta")
    
    ## call clustalw (needs to be installed and in the path!)
    system(paste(.findExecuable("clustalw -profile1="), prof1,
		    "-profile2=", prof2, " ", param2, " ", param, sep=""))

    read.DNAMultipleAlignment(paste(temp_file, ".aln", sep=""), 
	    format="clustal")
}
