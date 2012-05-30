\name{scoreSequence}
\alias{scoreSequence}
\alias{scoreSequence,GenModel,matrix-method}
\title{Score a New Sequence Given a Model}
\description{
Calculates a score of how likely it is that a new sequence was generated
by the same process as the sequences used to build the model.
}
\usage{
\S4method{scoreSequence}{GenModel,matrix}(x, newdata, 
    method = "supported_transitions", 
    match_cluster = "nn", plus_one = FALSE, 
    initial_transition = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ an \code{GenModel} object. }
  \item{newdata}{ a NSVSet containing the sequences to score. }
  \item{method}{ method to calculate the score (see details). }
  \item{match_cluster}{ do the new observations have to fall within
    the threshold of the cluster (\code{"exact"}) or is nearest neighbor 
    used (\code{"nn"})?}
  \item{plus_one}{ add one to each transition count. This is equal
      to start with one for each transition count, i.e. initially all
          transitions are equally likely. It prevents the product
          of probabilities to be zero if a transition was never observed. }
  \item{initial_transition}{ include the initial transition 
         in the computation?}
}
\details{
This function delegates the score computation to the implementation of
the model used. Currently only EMM models (see \pkg{rEMM}) are used.
The available scoring functions for these models are

\describe{
    \item{"product"}{
	Product of transition probabilities along the path of \eqn{x} in the 
	    model. A single missing transition (transition probability of zero) 
	    will result in
	    a score of 0. Use \code{plus_one} to avoid this.
	    \deqn{S_\mathrm{product} = \sqrt[l-1]{\prod_{i=1}^{l-1}{a_{s(i),s(i+1)}}}}{
		S_product = prod(a_s(i),s(i+1))^(1/(l-1))}
    }

    \item{"sum"}{
	Sum of transition probabilities along the path of \eqn{x} in the
	    model.
	    \deqn{S_\mathrm{sum} = \frac{1}{l-1} \sum_{i=1}^{l-1}{a_{s(i),s(i+1)}}}{
		S_sum = 1/(l-1) sum(a_s(i),s(i+1))}
    }
    \item{"log_sum"}{
	Sum of the log of the transition probabilities along the path of 
	    \eqn{x} in the model. The ranking of the scores is equivalent to
	    the product of probabilities, however, the calculation is more reliable 
	    since the product of probabilities might become a very small number.

	    A single missing transition (transition probability of zero)
	    will result in a score of neg. infinity. 
	    Use \code{plus_one} to avoid this. 

	    \deqn{S_\mathrm{log\_sum} = \frac{1}{l-1} \sum_{i=1}^{l-1}{\mathrm{log}(a_{s(i),s(i+1)})}}{
		S_sum = 1/(l-1) sum(log(a_s(i),s(i+1)))}
    }


    \item{"weighted_product"}{
	Weighted version of the product of probabilities. The weight is
	    the  similarity between a new data point and the state in the model
	    it is assigned to.
	    \deqn{S_\mathrm{weighted\_product} = \sqrt[l-1]{\prod_{i=1}^{l-1}{\mathrm{simil}(x_i,s(i))  a_{s(i),s(i+1)}}}}{
		P_weighted_product = prod(simil(x_i,s(i)) a_s(i),s(i+1))^(1/(l-1))}
    }
    \item{"weighted_sum"}{
	Weighted version of the sum of probabilities.
	    \deqn{S_\mathrm{weighted\_sum} = \frac{1}{l-1} \sum_{i=1}^{l-1}{\mathrm{simil}(x_i,s(i))  a_{s(i),s(i+1)}}}{
		S_weighted_sum = 1/(l-1) sum(simil(x_i,s(i))  a_s(i),s(i+1))}
    }
    \item{"weighted_log_sum"}{
	Weighted version of the sum of the log of probabilities.
	    \deqn{S_\mathrm{weighted\_log\_sum} = \frac{1}{l-1} \sum_{i=1}^{l-1}{\mathrm{log}(\mathrm{simil}(x_i,s(i))  a_{s(i),s(i+1)})}}{
		S_sum = 1/(l-1) sum(simil(x_i,s(i))  a_s(i),s(i+1))}
    }


    \item{"supported_transitions"}{Fraction of transitions in the new sequence \eqn{x} supported (present) in the model after assigning each data point in \eqn{x} to a state in
	the model.
	    \deqn{S_\mathrm{supported\_transitions} = \frac{1}{l-1} \sum_{i=1}^{l-1}{\mathrm{I}(a_{s(i),s(i+1)})}}{
		S_\mathrm{supported\_transitions} = 1/(l-1) sum(I(a_s(i),s(i+1)))}
    }
}

where 
\eqn{x_i} represents the \eqn{i}-th data point in the new sequence,
    \eqn{a(i,j)} is the transition probability from state \eqn{i}
    to state \eqn{j} in the model,
    \eqn{s(i)} is the state the \eqn{i}-th data point (\eqn{x_i}) in
    the new sequence is assigned to,
    and \eqn{\mathrm{simil}(\cdot)}{simil(.)} is a similarity function 
    (uses the same similarity/distance measure that was used to create the model;
     distances are converted into similarities using 
     \eqn{\mathrm{simil} = 1/(1+\mathrm{d})}{simil = 1/(1+d)}). 

    For missing transitions \eqn{\mathrm{I(v)}{I(v)}} is an indicator function which is 0 for \eqn{v=0} and 1 otherwise.
    }
\value{
    A vector with scores for the sequences.
}
\seealso{
    \code{\link{GenModel}} to create a new GenModel and
    \code{\link[rEMM]{score}} in \pkg{rEMM}.
 }
\examples{
db<-createGenDB("example.sqlite")
processSequences(system.file("examples/phylums", package="QuasiAlign"), db)

### create a model for Genus Syntrophomonas
model <- GenModelDB(db, rank="Genus", name="Syntrophomonas", table="NSV")
model

### score 5 random sequence against the model
random_sequence <- getSequences(db, table="NSV", limit=5, random=TRUE)
names(random_sequence)

scoreSequence(model, random_sequence)

scoreSequence(model, random_sequence, method ="product", plus_one=TRUE)

### clean up
closeGenDB(db)
unlink("example.sqlite")
}
\keyword{model}
