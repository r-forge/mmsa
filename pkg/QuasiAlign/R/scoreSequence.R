scoreSequence <- function(x, newdata, 
	method = "supported_transitions",
	match_cluster="weighted", prior=TRUE, ...) {


    if(is(newdata, "NSVSet")) {
	return(sapply(newdata, FUN=function(s) 
			scoreSequence(x, s, method, match_cluster, prior, ...)))
    }


    rEMM::score(x$model, newdata, method=method, 
                match_cluster=match_cluster, prior=prior, ...)
}

