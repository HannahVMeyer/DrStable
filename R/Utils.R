#' Print userinfo.
#'
#' Wrapper function around \code{\link{message}} that allows to turn the 
#' printing of messages to standard out.
#' on or off
#'
#' @param userinfo Vector of [string] element(s) and variables
#' @param verbose [boolean] If TRUE message is displayed on standard out, if 
#' FALSE, message is suppressed.
#' @param sep Delimiter [string] to separate message elements when userinfo 
#' given as vector.
#' @seealso \code{\link{message}} which this function wraps
#' @examples
#' # (not run)
#' # vmessage("Hello world", sep=' ')
vmessage <- function(userinfo, verbose=TRUE, sep=" ") {
    if (!is.character(sep)) {
        stop("Separator passed to vmessage must be of type character")
    }
    if (verbose) {
        message(paste(userinfo, collapse=sep))
    }
}

#' Sum of squares
#' 
#' Compute sum of squares of vector or matrix
#' 
#' @param data matrix of vector [float] to compute sum of squares for.
#' @return value [float] of sum of squares.
#' @examples 
#' x <- 1:10
#' x_ss <- ss(x)
ss <- function(data) {
    return(sum(data^2))
}

#' Inverse normalise
#' 
#' Apply rank-based quantile normalistation to data.
#' 
#' @param x vector with [n] elements to normalise.
#' @return x vector with [n] normalised elements.
#' @export
#' @examples 
#' x <- seq(1,20,2)
#' x_norm <- invnorm(x)
invnorm = function(x) {
    qnorm((rank(x, na.last="keep", ties.method="random") - 1.5)/sum(!is.na(x)))
}

#' Remove nulls from nested list
#' 
#' Remove NULL elements from list.
#' 
#' @param x [list] with [n] elements (NULL and non-NULL elements)
#' @return [list] with [n - \# NULL elements] elements
#' @export
#' @examples 
#' y <- list(a=1,b=2,c=NULL)
#' y_nnull <- rmnulls(y)
rmnulls <- function(x) {
    x[!vapply(x, is.null, logical(1))]
}

#' Comma-separated string to numeric vector.
#'
#' Split input of comma-separated string into vector of numeric, logical or
#' character values.
#'
#' @param commastring Input [character] vector containing numbers separated by
#' commas.
#' @param type Name [string] of the type of input variables; one of numeric
#' logical or character
#' @return Numeric vector of values extracted from commastring.
#' @export
commaList2vector <- function(commastring=NULL, type="numeric") {
    if (is.null(commastring)) {
        return(NULL)
    }
    commastring <- gsub(" ", "", commastring)
    tmp <- type.convert(unlist(strsplit(commastring, ",")), as.is=TRUE)
    if (any(!is.numeric(tmp)) && type == "numeric") {
        stop("Type numeric specified, but comma-list contains at least one
             non-numeric argument")
    }
    if (any(!is.logical(tmp)) && type == "logical") {
        stop("Type logical specified, but comma-list contains at least one
             non-logical argument")
    }
    if (any(!is.character(tmp)) && type == "character") {
        stop("Type character specified, but comma-list contains at least one
             non-character argument")
    }
    if (! type  %in% c("numeric", "logical", "character")) {
        stop("Unknown type of comma-separated list elements")
    }
    return(tmp)
}
