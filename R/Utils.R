vmessage <- function(userinfo, verbose=TRUE, sep=" ") {
    if (verbose) {
        message(paste(userinfo, collapse=sep))
    }
}

#' Sum of squares
ss <- function(data) {
    return(sum(data^2))
}

#' Inverse normalise
invnorm = function(x) {
    qnorm((rank(x, na.last="keep", ties.method="random") - 1.5)/sum(!is.na(x)))
}

#' Remove nulls from nested list
rmnulls <- function(x) {
    x[!vapply(x, is.null, logical(1))]
}