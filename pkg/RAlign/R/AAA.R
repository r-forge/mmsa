
.findExecuable <- function(exe) {
    path <- Sys.which(exe)
    if(path=="") stop("Executable for ", exe, " not found!", call.=FALSE)
    path
}

