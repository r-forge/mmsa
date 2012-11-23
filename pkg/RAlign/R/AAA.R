
.findExecuable <- function(exe) {
    path <- Sys.which(exe)
    if(all(path=="")) stop("Executable for ", paste(exe, collapse=" or "), " not found!", call.=FALSE)
    
    path[which(path!="")[1]]
}

