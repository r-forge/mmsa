
classifyRDP <- function(x, confidence=.8, property=NULL, java_args="-Xmx1g"){

    ## check 
    if(Sys.getenv("RDP_JAR_PATH") =="") stop("Environment variable 'RDP_JAR_PATH needs to be set!'")

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
    outfile <- paste(temp_file, "_tax_assignments.txt", sep="")
    
    ## property?
    if(!is.null(property)) property <- paste("-t", property)
    else property <- ""

    write.XStringSet(x, infile, append=FALSE, format="fasta")

    system(paste("java", java_args, "-jar", Sys.getenv("RDP_JAR_PATH"), 
		    "-f fixrank", property, "-q", infile, "-o", outfile),
	    ignore.stdout=TRUE, ignore.stderr=TRUE)

    ## read and parse rdp output
    cl_tab <- read.table(outfile, sep="\t") 
    
    ## remove empty columns
    cl_tab <- cl_tab[!sapply(cl_tab, FUN=function(x) all(is.na(x)))]
    
    seq_names <- cl_tab[,1] ## sequence names are in first column
    
    i <- seq(2, ncol(cl_tab), by=3) ## 3 columns for each tax. level
    
    ## get classification
    cl <- cl_tab[,i]	
    dimnames(cl) <- list(seq_names, as.matrix(cl_tab[1,i+1])[1,])

    ## get confidence
    conf <- as.matrix(cl_tab[,i+2])
    dimnames(conf) <- list(seq_names, as.matrix(cl_tab[1,i+1])[1,])
    
    if(confidence>0) cl[conf < confidence] <- NA
    
    attr(cl, "confidence") <- conf    
    cl
}

