#######################################################################
# RAlign - Interfaces to several sequence alignment and classification tools
# Copyright (C) 2012 Michael Hahsler and Anurag Nagar
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

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

    writeXStringSet(x, infile, append=FALSE, format="fasta")

    system(paste(.findExecuable("blastn"), "-db", db,
		    "-query", infile, "-out", outfile, "-outfmt 6", BLAST_args))

    ## read and parse rdp output
    cl_tab <- read.table(outfile)
    colnames(cl_tab) <- c( "QueryID",  "SubjectID", "Perc.Ident",
	    "Alignment.Length", "Mismatches", "Gap.Openings", "Q.start", "Q.end",
	    "S.start", "S.end", "E", "Bits" )
    
    cl_tab
}

