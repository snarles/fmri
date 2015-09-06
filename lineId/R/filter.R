#' Apply no filtering
#' 
#' A \code{filter_method}.
#' @export
no_filter <- function(X, Y, B, ...) rep(TRUE, dim(Y)[2])