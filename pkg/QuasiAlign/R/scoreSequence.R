scoreSequence <- function(x, newdata, 
	method = "supported_transitions",
	match_cluster="nn", plus_one = FALSE,
	initial_transition = FALSE) {


    if(is(newdata, "NSVSet")) {
	return(sapply(newdata, FUN=function(s) 
			scoreSequence(x, s, method, match_cluster, 
				plus_one, initial_transition)))
    }

    rEMM::score(x$model, newdata, method, match_cluster, 
	    plus_one, initial_transition)
}

