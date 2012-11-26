
BLAST <- function(x, db=NULL, BLAST_args="") {

    if(is.null(db)) stop("Specify the database as argument db!")

    ## get temp files and change working directory
    wd <- tempdir()
    dir <- getwd()
    temp_file <- basename(tempfile(tmpdir = wd))
    on.exit({
		file.remove(Sys.glob(paste(temp_file, "*", sep="")))
		setwd(dir)
	    })
    setwd(wd)

    infile <- paste(temp_file, ".fasta", sep="")
    outfile <- paste(temp_file, "_BLAST_out.txt", sep="")

    write.XStringSet(x, infile, append=FALSE, format="fasta")

    system(paste(.findExecuable("blastn"), "-db", db,
		    "-query", infile, "-out", outfile, "-outfmt 6", BLAST_args))

    ## read and parse rdp output
    cl_tab <- read.table(outfile)
    colnames(cl_tab) <- c( "QueryID",  "SubjectID", "Perc.Ident",
	    "Alignment.Length", "Mismatches", "Gap.Openings", "Q.start", "Q.end",
	    "S.start", "S.end", "E", "Bits" )
    
    cl_tab
}

