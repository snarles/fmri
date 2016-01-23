#' Sum of correct posterior probabilities
#' 
#' A \code{scoring_method}.
#' @param plikes Matrix of posterior likelihoods, rows = responses, columns = classes
#' @param i_chosen True classes
#' @return Sum of correct probabilities
#' @export
sum_score <- function(plikes, i_chosen, ...) {
  mat <- apply(plikes, 1,  function(v) {
    v <- v - max(v)
    v <- exp(v)
    v/sum(v)
  })
  sum(t(mat)[cbind(1:dim(mat)[1], i_chosen)])
}

#' Counting correct classifications within top-K rankings
#' 
#' A \code{scoring_method}.
#' @param plikes Matrix of posterior likelihoods, rows = responses, columns = classes
#' @param i_chosen True classes
#' @param k Number of top candidates for each response
#' @return Number correct
#' @export
topk_score <- function(plikes, i_chosen, k = 1, ...) {
  mmat <- cbind(i_chosen, plikes)
  res <- apply(mmat, 1, function(v) v[1] %in% order(-v[-1])[1:k])
  sum(res)
}

#' Probability of misclassification if we resampled m out of K classes
#' 
#' A \code{scoring_method}.
#' @param plikes Matrix of posterior likelihoods, rows = responses, columns = classes
#' @param i_chosen True classes
#' @param m number of classes resampled
#' @param replace TRUE or FALSE: whether to sample with replacement
#' @return Misclassification probability
#' @export
resample_misclassification <- function(plikes, i_chosen,
                                       m = dim(plikes)[2], 
                                       replace = FALSE, ...)
{
  K <- dim(plikes)[2]
  mmat <- cbind(i_chosen, plikes)
  nbads <- apply(mmat, 1, function(v) sum(v[-1] >= v[-1][v[1]]))
  p_i <- 1 - nbads/K
  if (replace) {
    ans <- 1 - mean(p_i ^ (m - 1))    
  } else {
    ans <- 1 - mean(dhyper(0, nbads - 1, K - nbads, m - 1))
  }
  ans
}

