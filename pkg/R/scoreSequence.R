scoreSequence <- function(x, newdata, 
	method = "product",
	match_cluster="nn", plus_one = FALSE,
	initial_transition = FALSE)
{
	rEMM::score(x$model, newdata, method, match_cluster, plus_one, initial_transition)
}

