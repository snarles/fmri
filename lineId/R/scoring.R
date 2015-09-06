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