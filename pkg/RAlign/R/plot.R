### plot as a Sequence Logo

### DNA
plot.DNAMultipleAlignment <- function(x, start=1, end=NULL, ic.scale=FALSE, ...) plot(as(x, "DNAStringSet"), start, end, ic.scale, ...)


plot.DNAStringSet <- function(x, start=1, end=NULL, ic.scale=FALSE, ...)
{
    cm <- consensusMatrix(x,baseOnly=TRUE)
    if (is.null(end))
	end <- ncol(cm)
    cm <- cm[1:4,start:end]
    for(i in 1:ncol(cm))
	cm[,i] <- cm[,i]/colSums(cm)[i]
    pwm <- seqLogo::makePWM(cm)
    seqLogo::seqLogo(pwm, ic.scale,  ...)
}

### RNA
### FIXME: SeqLogo only supports DNA
plot.RNAMultipleAlignment <- function(x, start=1, end=NULL, ic.scale=FALSE, ...) plot(as(x, "RNAStringSet"), start, end, ic.scale, ...)


plot.RNAStringSet <- function(x, start=1, end=NULL, ic.scale=FALSE, ...)
{
    warning("Sequence logo only supports DNA. All U's are plotted as T's!")
    cm <- consensusMatrix(x,baseOnly=TRUE)
    if (is.null(end))
	end <- ncol(cm)
    cm <- cm[1:4,start:end]
    for(i in 1:ncol(cm))
	cm[,i] <- cm[,i]/colSums(cm)[i]
    pwm <- seqLogo::makePWM(cm)
    seqLogo::seqLogo(pwm, ic.scale,  ...)
}

