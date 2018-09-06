#' Command line execution for Dimensionality reduction
#' 
#' runCreateSubsets runs without arguments. Upon call, it reads 
#' command-line parameters and supplies these to \code{\link{}} and 
#' \code{\link{}}. For details on input to \code{\link{}} 
#' and \code{\link{}}, please refer to their help pages. For help on 
#' the command line arguments that can be passed, see examples below. From the 
#' command line, the help function can be called via 
#' 'Rscript -e "Stability::runCreateSubsets()" --args --help 
#'  
#' @export
#' 
#' @examples
#' # (not run)
#' # :
#' # (not run)
#' # Rscript -e "Stability::runCreateSubsets()" \
#' #--args \ 
#' #--showProgress \
runCreateSubsets <- function() {
    option_list <- list(
        make_option(c("-hd", "--highdim"), action="store", dest="highdimfile", 
                    type="character", help="Path to comma-separated input file 
                    with N samples (rows) and P traits (columns). Subsetting 
                    occurs along rows, creating [M x P] subsets of the original
                    [N x P] dataset [default: %default]."),
        make_option(c("-c", "--cov"), action="store", dest="covariatefile", 
                    type="character", default=NULL, 
                    help="Path to (optional) comma-separated covariate file 
                    with N samples (rows) and C covariates (columns). [default: 
                    %default]."),
        make_option(c("--seed"), action="store", dest="seed", default=234,
                    type="integer", help="Seed to initialise random number 
                    generator [default: %default]."),
        make_option(c("-n", "--nrCV"), action="store", dest="nrCV", 
                    default=10, type="integer", help="Number of cross validation
                    sets to generate [default: %default]."),
        make_option(c("-s", "--sizeCV"), action="store", dest="sizeCV", 
                    default=0.8, type="double", help="Proportion of sampples in 
                    cross validation data sets [default: %default]."),
        make_option(c("--showProgress"), action="store_true", dest="verbose", 
                    default=FALSE, type="logical", help="If set, progress 
                    messages about simulation steps are printed to standard out
                    [default: %default]."))
    
    if (args$verbose) message("Parse command line arguments")
    args <- parse_args(OptionParser(option_list=option_list))
    
    # phenotype: N x P matrix
    if (!grepl("csv", args$highdimfile)) {
        stop("High-dimensional data files is not .csv")
    }
    if (args$verbose) {
        message("Read high-dimensional dataset from ", args$highdimfile)
    }
    Y <-  data.table::fread(args$highdimfile, sep=",", header=TRUE,
                            stringsAsFactors=FALSE, data.table=FALSE)
    rownames(Y) <- Y[,1]
    Y <- Y[,-1]
    
    if(!is.numeric(Y)) {
        stop("High-dimensional dataset contains at least one non-numeric entry")
    }
    
    # potential covariates
    if (!is.null(args$covariatefile)) {
        # N x C matrix
        if (!grepl("csv", args$covariatefile)) {
            stop("Covariates files is not .csv")
        }
        if (args$verbose) {
            message("Read covariates from ", args$covariatefile)
        }
        covs <- data.table::fread(args$covariatefile, sep=",", header=TRUE,
                                  stringsAsFactors=FALSE, data.table=FALSE)
        if(!is.numeric(covs)) {
            stop("Covarariates contain at least one non-numeric entry")
        }
        if (args$verbose) message("Regress covariates...")
        Y <- lm(Y ~ covs)$residuals
    }
    
    if (args$verbose) message("Set seed to ", seed)
    set.seed(seed)
    sample_matrix <- sapply(1:args$nrCV, function(x) {
        sample(nrow(Y), args$sizeCV * nrow(Y))
        })
    if (args$verbose) {
        message("Save cross-validation datasets to ", 
                gsub(".csv", "", args$highdimfile), "_cv...")
    dr <- lapply(1:ncol(sample_matrix), function(x){
            y_cv <- Y[sample_matrix[,x],]
            cvfile <- paste(gsub(".csv", "", args$highdimfile), "_cv", x, ".csv",
                            sep="")
            write.csv(y_cv, cvfile, quote=FALSE, row.names=TRUE, col.names=NA,
                      sep=",")
    })
}

#' Command line execution for Dimensionality reduction
#' 
#' runDimensionalityReduction runs without arguments. Upon call, it reads 
#' command-line parameters and supplies these to \code{\link{}} and 
#' \code{\link{}}. For details on input to \code{\link{}} 
#' and \code{\link{}}, please refer to their help pages. For help on 
#' the command line arguments that can be passed, see examples below. From the 
#' command line, the help function can be called via 
#' 'Rscript -e "Stability::runDimensionalityReduction()" --args --help 
#'  
#' @export
#' 
#' @examples
#' # (not run)
#' # :
#' # (not run)
#' # Rscript -e "Stability::runDimensionalityReduction()" \
#' #--args \ 
#' #--directory=/tmp \
#' #--showProgress \
runDimensionalityReduction <- function() {
    option_list <- list(
        make_option(c("-hd", "--highdim"), action="store", dest="highdimfile", 
                    type="character", help="Path to comma-separated input file 
                    with N samples (rows) and P traits (columns). Subsetting 
                    occurs along rows, creating [M x P] subsets of the original
                    [N x P] dataset [default: %default]."),
        make_option(c("-c", "--cov"), action="store", dest="covariatefile", 
                    type="character", default=NULL, 
                    help="Path to (optional) comma-separated covariate file 
                    with N samples (rows) and C covariates (columns). [default: 
                    %default]."),
        make_option(c("-s", "--seed"), action="store", dest="seed", default=234,
                    type="integer", help="Seed to initialise random number 
                    generator [default: %default]."),
        make_option(c("-n", "--nrCV"), action="store", dest="nrCV", 
                    default=10, type="integer", help="Number of cross validation
                    sets to generate [default: %default]."),
        make_option(c("-s", "--sizeCV"), action="store", dest="sizeCV", 
                    default=0.8, type="double", help="Proportion of samples in 
                    cross validation data sets [default: %default]."),
        make_option(c("-dim", "--dimensions"), action="store", dest="dim", 
                    default=NULL, type="integer", help="Maximum dimensionality 
                    [int] to retain in the data; large values can cause long
                    computation times; if not provided max(P,N) is chosen 
                    [default: %default]."),
        make_option(c("--neighbours"), action="store", dest="nOpt", 
                    default=NULL, type="integer", help="Number of neighbours
                    considered for dimensionality reduction; input for LLE,
                    LLE, LaplacianEigenmaps, Isomap, tSNE; if not provided, 
                    will be estimated via lle::calc_k. For details see 
                    'computeDimReduction' function [default: %default]."),
        make_option(c("--kmax"), action="store", dest="dim", 
                    default=40, type="integer", help="if neighbours is not 
                    provided, kmax [int] specifies the maximum number of 
                    neighbours supplied to lle::calc_k [default: %default]."),
        make_option(c("--kmin", "--dimensions"), action="store", dest="dim", 
                    default=1, type="integer", help="if neighbours is not 
                    provided, kmin [int] specifies the minimum number of 
                    neighbours supplied to lle::calc_k [default: %default]."),
        make_option(c("-m", "--method"), action="store", dest="method", 
                    type="character", default=NULL, 
                    help="Dimensionality reduction method to apply 
                    [default: %default]."),
        make_option(c("-o", "--outdir"), action="store", dest="outdir", 
                    type="character", default=NULL, 
                    help="Path to directory where dimensionality reduction
                    results will be saved [default: %default]."),
        make_option(c("--showProgress"), action="store_true", dest="verbose", 
                    default=FALSE, type="logical", help="If set, progress 
                    messages about simulation steps are printed to standard out
                    [default: %default]."))
    
    if (args$verbose) message("Parse command line arguments")
    args <- parse_args(OptionParser(option_list=option_list))
    
    # phenotype: N x P matrix
    if (!grepl("csv", args$highdimfile)) {
        stop("High-dimensional data files is not .csv")
    }
    if (args$verbose) {
        message("Read high-dimensional dataset from ", args$highdimfile)
    }
    Y <-  data.table::fread(args$highdimfile, sep=",", header=TRUE,
                            stringsAsFactors=FALSE, data.table=FALSE)
    rownames(Y) <- Y[,1]
    Y <- Y[,-1]
    
    if(!is.numeric(Y)) {
        stop("High-dimensional dataset contains at least one non-numeric entry")
    }
    
    # potential covariates
    if (!is.null(args$covariatefile)) {
        # N x C matrix
        if (!grepl("csv", args$covariatefile)) {
            stop("Covariates files is not .csv")
        }
        if (args$verbose) {
            message("Read covariates from ", args$covariatefile)
        }
        covs <- data.table::fread(args$covariatefile, sep=",", header=TRUE,
                                  stringsAsFactors=FALSE, data.table=FALSE)
        if(!is.numeric(covs)) {
            stop("Covarariates contain at least one non-numeric entry")
        }
        if (args$verbose) message("Regress covariates...")
        Y <- lm(Y ~ covs)$residuals
    }
    
    if (is.null(args$dim)) {
        
    }
    if (crossvalidate) {
        if (args$verbose) message("Set seed to ", seed)
        set.seed(seed)
        sample_matrix <- sapply(1:args$nrCV, function(x) {
            sample(nrow(Y), args$sizeCV * nrow(Y))
        })
        if (args$verbose) {
            message("Save cross-validation datasets to ", 
                    gsub(".csv", "", args$highdimfile), "_cv...")
        }
        dr <- lapply(1:ncol(sample_matrix), function(x){
            y_cv <- Y[sample_matrix[,x],]
            cvfile <- paste(gsub(".csv", "", args$highdimfile), "_cv", x,
                            ".csv", sep="")
            write.csv(y_cv, cvfile, quote=FALSE, row.names=TRUE, col.names=NA,
                          sep=",")
            output_cv <- paste(args$outdir, "/CV", x, "_", sep="")
            if (verbose) {
                message("Dimensionality reduction for crossvalidation:", x)
            }
            dr <- computeDimReduction(m=args$method, Y=y_cv, output=output_cv,
                                      ndim=args$dim, optN=args$optN,
                                      kmin=args$kmin, kmax=args$kmax)
        })
    } else {
        if (gsub(".*(.)", "\\1", output) != "_") {
            output <- paste(output, "/", sep="")
        }
        dr <- computeDimReduction(m=args$method, Y=Y, output=output,
                                  ndim=args$dim)
    }
}