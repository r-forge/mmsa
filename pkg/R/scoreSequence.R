scoreSequence <- function(x, newdata, 
	method = c("prod", "malik", "misstran","compareTransitions"),
	match_cluster="nn", plus_one = FALSE,
	initial_transition = FALSE)
{
            #for Kullback method add 1 to newdata
            if(x$model@measure=="Kullback") newdata <- newdata+1

            method <- match.arg(method)

            if(method == "prod") {
				sc<-rEMM::score(x$model, newdata=newdata, method="prod",
                    match_cluster=match_cluster, plus_one=plus_one,
                    initial_transition=initial_transition)
				return(sc)
            }

            if(method == "malik") {
                emm_test<- EMM(threshold=attr(x$model,"threshold"))
                emm_test<-build(emm_test,newdata)
                distance<-0
                d<-dist(cluster_centers(emm_test),cluster_centers(x$model))
                closestClusters<-apply(d, MARGIN=1, FUN=which.min)
                for(i in 1:length(closestClusters))
                {
                        closestDist <- d[i,closestClusters[i]]
                        #check if previous transition exists in model
                        previousTransition<-transition(x$model,as.character(closestClusters[i]),as.character(closestClusters[i]-1))
                        if (previousTransition > 0)
                            previousTransition=1
                        else if (previousTransition ==0)
                            previousTransition = -log(1/nstates(emm_test));

                        distance<-distance+(closestDist*previousTransition)
                }
            return (1/(1+distance))

            }

            if(method == "misstran") {
                transitionTable <- transition_table(x$model, newdata, method="prob",
                        match_cluster, plus_one=FALSE,
                        initial_transition)
                missingTransitions <- sum(transitionTable[,3]==0)
                return(missingTransitions)
            }

            if(method == "compareTransitions") {
                transitionTable <- transition_table(x$model, newdata, method="prob",
                        match_cluster, plus_one=TRUE,
                        initial_transition)
				dist<-0
				for(i in 1:nrow(transitionTable))
				{	
					tempDist<-dist(cluster_centers(x)[transitionTable[i,1]],cluster_centers(x)[transitionTable[i,2]],measure=x$measure)
					tempSim <- pr_dist2simil(tempDist)
					tempProb <- transitionTable[i,3]
					dist<-dist+(abs(log(tempSim))*tempProb)
				}
				return(1/(1+dist))
            }


}

