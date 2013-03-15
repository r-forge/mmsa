
.onAttach <- function(lib, pkg) {
  registerDoSEQ()
}
  
.findExecuable <- function(exe) {
    path <- Sys.which(exe)
    if(all(path=="")) stop("Executable for ", paste(exe, collapse=", "), " not found!", call.=FALSE)
    
    path[which(path!="")[1]]
}
