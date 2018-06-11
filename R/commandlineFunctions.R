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
    cat(unlist(args), "\n")
    
    # path to input file (formatted via Phenotype_setup.R)
    inputfile <- args$inputfile
    # path to covariate file
    covariatefile <- args$covariatefile 
    # output_prefix
    output <- args$output
    # method(s)
    method <-  unlist(strsplit(args$method, ","))
    # cv
    crossvalidate <- args$cv
    # seed
    seed <- as.numeric(args$seed)
    # seed
    ndim <- as.numeric(args$ndim)
    # only run to subset data
    createSubsets <- as.logical(args$createSubsets)
    
    ################
    ### analysis ###
    ################
    
    # phenotype: N x P matrix
    Y <-  readRDS(inputfile)
    
    # potential covariates
    if (!is.null(covariatefile)) {
        # N x C matrix
        covs <- readRDS(covariatefile)
        Y <- lm(Y ~ covs)$residuals
    }
    
    if (crossvalidate) {
        set.seed(seed)
        sample_matrix <- sapply(1:10, function(x) sample(nrow(Y), 4/5*nrow(Y)))
        dr <- lapply(1:ncol(sample_matrix), function(x){
            y_cv <- Y[sample_matrix[,x],]
            output_cv <- paste(output, "/CV", x, "_", sep="")
            cvfile <- paste(gsub(".rds", "", inputfile), "_cv", x, ".rds", sep="")
            saveRDS(y_cv, cvfile)
            if (! createSubsets) {
                vmessage(c("Crossvalidation:", x))
                dr <- computeDimReduction(m=method, Y=y_cv, output=output_cv, 
                                          ndim=ndim)
            }
        })
    } else {
        if (gsub(".*(.)", "\\1", output) != "_") {
            output <- paste(output, "/", sep="")
        }
        dr <- computeDimReduction(m=method, Y=Y, output=output, ndim=ndim)
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
#' 'Rscript -e "Stability::estimateStability()" --args --help 
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