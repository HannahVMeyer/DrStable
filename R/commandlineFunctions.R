#' Command line execution for Dimensionality reduction
#'
#' runCreateSubsets runs without arguments. Upon call, it reads
#' command-line parameters and supplies these to \code{\link{}} and
#' \code{\link{}}. For details on input to \code{\link{}}
#' and \code{\link{}}, please refer to their help pages. For help on
#' the command line arguments that can be passed, see examples below. From the
#' command line, the help function can be called via
#' Rscript -e "Stability::runCreateSubsets()" --args --help
#'
#' @export
#'
#' @examples
#' # (not run)
#' # :
#' # (not run)
#' # Rscript -e "Stability::runCreateSubsets()" \
#' #--args
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

    args <- parse_args(OptionParser(option_list=option_list))
    if (args$verbose) message("Parse command line arguments")

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
    Y <- as.matrix(Y[,-1])

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
        rownames(covs) <- covs[,1]
        covs <- as.matrix(covs[,-1])

        if(!is.numeric(covs)) {
            stop("Covariates contain at least one non-numeric entry")
        }
        if (!all(rownames(covs) == rownames(Y))) {
            stop("Rownames of highdim file and covariate file differ")
        }
        if (args$verbose) message("Regress covariates...")
        Y <- lm(Y ~ covs)$residuals
    }

    if (args$verbose) message("Set seed to ", args$seed)
    set.seed(args$seed)
    sample_matrix <- sapply(1:args$nrCV, function(x) {
        sample(nrow(Y), args$sizeCV * nrow(Y))
        })
    if (args$verbose) {
        message("Save cross-validation datasets to ",
                gsub(".csv", "", args$highdimfile), "_cv...")
    }
    dr <- lapply(1:ncol(sample_matrix), function(x){
            y_cv <- Y[sample_matrix[,x],]
            cvfile <- paste(gsub(".csv", "", args$highdimfile), "_cv", x, ".csv",
                            sep="")
            write.table(y_cv, cvfile, quote=FALSE, row.names=TRUE, col.names=NA,
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
#' Rscript -e "Stability::runDimensionalityReduction()" --args --help
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
        make_option(c("--seed"), action="store", dest="seed", default=234,
                    type="integer", help="Seed to initialise random number
                    generator [default: %default]."),
        make_option(c("-cv", "--crossvalidate"), action="store_true",
                    dest="crossvalidate", default=FALSE, help="Split dataset in
                    nrCV of size sizeCV and apply dimreduction to these subsets
                    [default: %default]."),
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
        make_option(c("--neighbours"), action="store", dest="optN",
                    default=NULL, type="integer", help="Number of neighbours
                    considered for dimensionality reduction; input for LLE,
                    LLE, LaplacianEigenmaps, Isomap, tSNE; if not provided,
                    will be estimated via lle::calc_k. For details see
                    'computeDimReduction' function [default: %default]."),
        make_option(c("--kmax"), action="store", dest="kmax",
                    default=40, type="integer", help="if neighbours is not
                    provided, kmax [int] specifies the maximum number of
                    neighbours supplied to lle::calc_k [default: %default]."),
        make_option(c("--kmin"), action="store", dest="kmin",
                    default=1, type="integer", help="if neighbours is not
                    provided, kmin [int] specifies the minimum number of
                    neighbours supplied to lle::calc_k [default: %default]."),
        make_option(c("-m", "--method"), action="store", dest="method",
                    type="character", default=NULL,
                    help="Dimensionality reduction method to apply
                    [default: %default]."),
        make_option(c("--name"), action="store", dest="name",
                    type="character", default=NULL, help="Name of output file
                    [default: %default]."),
        make_option(c("-o", "--outdir"), action="store", dest="outdir",
                    type="character", default=NULL,
                    help="Full path/to/directory [string] where results will be
                    saved. Full file name will be created by
                    outdir/name_method.csv. Alternatively, prefix can be
                    provided to create file prefix_method.csv [default:
                    %default]."),
        make_option(c("-p", "--prefix"), action="store", dest="prefix",
                    type="character", default=NULL, help="Prefix [string] of
                    output file name, containing full path/to/directory [string]
                    where results will be saved. Results will be saved
                    as prefix_method.csv  %default]."),
        make_option(c("--showProgress"), action="store_true", dest="verbose",
                    default=FALSE, type="logical", help="If set, progress
                    messages about simulation steps are printed to standard out
                    [default: %default]."))

    args <- parse_args(OptionParser(option_list=option_list))
    if (args$verbose) message("Parse command line arguments")

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
    Y <- as.matrix(Y[,-1])

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
        rownames(covs) <- covs[,1]
        covs <- as.matrix(covs[,-1])

        if(!is.numeric(covs)) {
            stop("Covarariates contain at least one non-numeric entry")
        }
        if (!all(rownames(covs) == rownames(Y))) {
            stop("Rownames of highdim file and covariate file differ")
        }
        if (args$verbose) message("Regress covariates...")
        Y <- lm(Y ~ covs)$residuals
    }

    if (args$crossvalidate) {
        if (args$verbose) message("Set seed to ", args$seed)
        set.seed(args$seed)
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
            write.table(y_cv, cvfile, quote=FALSE, row.names=TRUE, col.names=NA,
                          sep=",")
            output_cv <- paste(args$outdir, "/CV", x, "_", sep="")
            if (verbose) {
                message("Dimensionality reduction for crossvalidation:", x)
            }
            dr <- computeDimReduction(m=args$method, Y=y_cv,
                                      verbose=args$verbose,
                                      ndim=args$dim, optN=args$optN,
                                      kmin=args$kmin, kmax=args$kmax)
            saveDimReduction(results=dr, prefix=args$prefix, outdir=args$outdir,
                             name=args$name, method=args$method)
        })
    } else {
        dr <- computeDimReduction(m=args$method, Y=Y,
                                      verbose=args$verbose,
                                      ndim=args$dim, optN=args$optN,
                                      kmin=args$kmin, kmax=args$kmax)
        saveDimReduction(results=dr, prefix=args$prefix, outdir=args$outdir,
                         name=args$name, method=args$method)
    }
}


#' Command line execution for Stability
#'
#' estimateStability runs without arguments. Upon call, it reads command-line
#' parameters and supplies these to \code{\link{}} and
#' \code{\link{}}. For details on input to \code{\link{}}
#' and \code{\link{}}, please refer to their help pages. For help on
#' the command line arguments that can be passed, see examples below. From the
#' command line, the help function can be called via
#' Rscript -e "Stability::estimateStability()" --args --help
#'
#' @export
#'
#' @examples
#' # (not run)
#' #
#' # (not run)
#' # Rscript -e "Stability::estimateStability()" \
#' #--args \
#' #--directory=/tmp \
#' #--showProgress \


runEstimateStability <- function() {
    option_list <- list(
        make_option(c("--NrSamples"), action="store", dest="NrSamples",
                    type="integer", help="Number of samples to
                    simulate [default: %default]."),
        make_option(c("--NrPhenotypes"), action="store", dest="NrPhenotypes",
                    type="integer", help="Number of phenotypes to
                    simulate [default: %default]."),
        make_option(c("--showProgress"), action="store_true", dest="verbose",
                    default=FALSE, type="logical", help="If set, progress
                    messages about simulation steps are printed to standard out
                    [default: %default]."))

    args <- parse_args(OptionParser(option_list=option_list))
    args <- commandArgs(asValue= TRUE, defaults=list(cv=FALSE, seed=234,
                                                     ndim=NULL, covariatefile=NULL,
                                                     createSubsets=FALSE))
    args <- commandArgs(asValue= TRUE,
                        defaults=list(cv=10))

    cat(unlist(args), "\n")

    # path to input file
    analysisdir <- args$dir
    # NrSamples
    N <- as.numeric(args$N)
    # output_prefix
    output <- args$output
    # method(s)
    method <-  args$method
    # cv
    crossvalidate <- as.numeric(args$cv)

    dr  <- lapply(1:crossvalidate, function(cv, d) {
        filecv <- paste(d, "/CV", cv, "_",  method, "_reducedY.csv",
                        sep="")
        if (file.exists(filecv)) {
            datacv <- t(as.matrix(read.csv(filecv, sep=",", header=TRUE,
                                           row.names=1)))
            if (method == "") {
                datacv <- datacv[-1,]
                rownames(datacv) <- paste("DR", 1:100, sep="")
            }
        } else {
            datacv <- NULL
        }
        return(datacv)
    }, d=analysisdir)
    names(dr) <- paste("cv", 1:crossvalidate, sep="")

    dr <- rmnulls(dr)

    if (length(dr) > 1) {
        cat("Determine robustness of ", method, " for ", N, " samples")
        color=c(wes_palette(20, name="Zissou", type='continuous'))
        threshold=seq(0.0, 0.95, 0.05)
        results <- analyseDimreduction(dr, output=output, threshold=threshold)
        saveRDS(results, paste(output, "/", method, "_robustness.rds", sep=""))
    } else {
        cat("Not enough datasets (#", length(dr), ") to compute robustness")
    }
}
