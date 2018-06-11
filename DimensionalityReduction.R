#################
### libraries ###
#################

library ("R.methodsS3")
library ("R.oo")
library ("R.utils")
library("methods")

#################
### functions ###
#################

source("~/GWAS/analysis/dimred/Trustworthiness.R")

vmessage <- function(userinfo, verbose=TRUE, sep=" ") {
    if (verbose) {
        message(paste(userinfo, collapse=sep))
    }
}

dimReduction <- function(Y, ndim, method=c("DiffusionsMaps", "DRR", "ICA", 
                                           "LLE", "Isomap", "LaplacianEigenmap", 
                                            "MDS", "PCA","kPCA", "nMDS", 
											"tSNE", "PEER"),
                         optN=NULL) {
    if (method == "DiffusionsMaps") {
        results <- diffusionMap::diffuse(dist(Y), neigen=ndim)
        reducedY <- results$X
    } else if (method == "DRR") {
        results <- DRR::drr(Y, ndim=ndim)
        reducedY <- results$fitted.data
    }  else if (method == "ICA") {
        results <- fastICA::fastICA(Y, n.comp=ndim, fun= "logcosh", method="C")
        #results <- ica::icajade(Y, nc=ndim)
        reducedY <- results$S
    } else if (method == "LLE") {
        results <- lle::lle(Y, k=optN, m=ndim, id=TRUE)
        reducedY <- results$Y
    } else if (method == "Isomap") {
        results <- vegan::isomap(dist(Y), ndim=100, k=optN, fragmentedOK=TRUE)
        reducedY <- results$points
    } else if (method == "LaplacianEigenmap") {
        library("dimRed")
        y <- as(as.data.frame(Y), "dimRedData")
        results <- dimRed::embed(y, 'LaplacianEigenmaps', ndim=100, knn=optN)
        reducedY <- results@data@data
    } else if (method == "MDS") {
        results <- stats::cmdscale(dist(Y), k=ndim)
        reducedY <- results  
    } else if (method == "PCA") {
        results <- stats::prcomp(Y)
        reducedY <- results$x[,1:ndim]
    } else if (method == "kPCA") {
        results <- kernlab::kpca(Y, features=ndim)
        reducedY <- results@rotated
    } else if (method == "nMDS") {
        results <- vegan::metaMDS(dist(Y), k=ndim)
        reducedY <- results$points
	} else if (method == "tSNE") {
        results <- Rtsne::Rtsne(Y, dims=ndim, initial_dims=ndim, 
                                perplexity=optN)
        reducedY <- results$Y
	} else if (method == "PEER") {
    	model = peer::PEER()
    	# Set observed data
    	peer::PEER_setPhenoMean(model, Y)
        peer::PEER_setAdd_mean(model, TRUE)
    	peer::PEER_setNk(model, 100)
    	peer::PEER_setNmax_iterations(model, 1000)
    	peer::PEER_update(model)
    	results <- list(Y=peer::PEER_getX(model)[,-1], W=peer::PEER_getW(model), 
					precision=peer::PEER_getAlpha(model))
		reducedY <- results$Y
    } else {
        stop("Method: ", method, " does not exist, possible methods are: ",
             "DiffusionsMaps, DRR, ICA, LLE, Isomap, LaplacianEigenmap, MDS,",
			  "PCA,", "kPCA, nMDS, tSNE and PEER")
    }
    return(list(results=results, reducedY=reducedY))
}

computeDimReduction <-  function(method, Y, output, optN=NULL, ndim=NULL) {
    # phenotype dimensions
    N <- nrow(Y)
    P <- ncol(Y)
    
    # number of dimensions to estimate
    if (is.null(ndim)) {
        if (P < N) {
            ndim <- P
        } else {
            ndim <- N
        }
    }
    if (any(method %in% c("LLE", "LaplacianEigenmaps", "Isomap", "tSNE"))) {
        nEstimate <- ndim/10
        if (nEstimate > 100) nEstimate = 100 
        library("lle")
        # find neighbours
        neighbours <- calc_k(Y, m=nEstimate, kmax=40, parallel=TRUE,
                                                         plotres=FALSE)
        optN <- neighbours$k[which.min(neighbours$rho)]
    }
    lapply(method, function(m, Y, output, optN, ndim) {
        vmessage(c("Running dimensionality reduction:", m))
        # dimensionalityReduction
        red <- dimReduction(Y=Y, method=m, ndim=ndim, optN=optN)
        rownames(red$reducedY) <- rownames(Y)
        colnames(red$reducedY) <- paste("DR", 1:ncol(red$reducedY), sep="")
        
        neighbours <- round(1:5 * nrow(Y)/100)
        M <- TandC(original=Y, reduced=red$reducedY, 
                   neighbours=neighbours)
        
        saveRDS(red$results, paste(output,  m, "_results.rds", sep=""))
        write.table(red$reducedY, paste(output,  m, "_reducedY.csv", 
                                        sep=""), 
                    quote=FALSE, row.names=TRUE, col.names=NA, sep=",")
        write.table(M, paste(output, m, "_TandC.csv", 
                                        sep=""), 
                    quote=FALSE, row.names=TRUE, col.names=NA, sep=",")
        return(M)
    }, Y=Y, output=output, optN=optN, ndim=ndim)
}


############
### data ###
############

### read command line arguments
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




