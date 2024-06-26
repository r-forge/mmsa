#######################################################################
# BiostringsTools - Interfaces to several sequence alignment and 
# classification tools
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

blast  <- function(db = NULL) {
  if(is.null(db)) stop("No BLAST database specified!")
  db <- file.path(normalizePath(dirname(db)), basename(db))
  if(length(Sys.glob(paste(db, "*", sep="")))<1) stop("BLAST database does not exit!")
  
  structure(list(db = db), class="BLAST")
}

print.BLAST <- function(x, info=TRUE, ...) {
  cat("BLAST Database\nLocation:", x$db, "\n")
  
  if(info) {
    out <- system(paste(.findExecutable("blastdbcmd"), "-db", x$db,
      "-info"), intern=TRUE)
    cat(paste(out, collapse="\n"))
    cat("\n")
  }
}

blast_help <- function() {
  system(paste(.findExecutable(c("blastn")),
    "-help"))
}


predict.BLAST <- function(object, newdata, BLAST_args="", ...) {
  
  db <- object$db
  x <- newdata
  
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
  
  system(paste(.findExecutable("blastn"), "-db", db,
    "-query", infile, "-out", outfile, "-outfmt 6", BLAST_args))
  
  ## read and parse rdp output
  cl_tab <- read.table(outfile)
  colnames(cl_tab) <- c( "QueryID",  "SubjectID", "Perc.Ident",
    "Alignment.Length", "Mismatches", "Gap.Openings", "Q.start", "Q.end",
    "S.start", "S.end", "E", "Bits" )
  
  cl_tab
}



