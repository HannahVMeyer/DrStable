#' Save results of dimensionality reduction
#'
#' @param results named list of results from dimensionality reduction;
#' Y_red: named list with dimensionality reduced phenotypes (reducedY) and
#' object returned by specified dimensionality reduction method (results) with
#' additional output.
#' M: vector [double] with Trustworthiness and Continuity estimates for the
#' dimensionality reduction.
#' @param name [string] name of output file
#' @param outdir [string] full path/to/directory [string] where results will be
#' saved. Full file name will be created by outdir/name_method.csv.
#' Alternatively, prefix can be provided to create file prefix_method.csv.
#' @param prefix [string] of output file name, containing full
#' path/to/directory [string] where results will be saved. Results will be saved
#' as prefix_method.csv
#' @param  method [string] name of dimensionality reduction method; will be
#' appended to output file name.
#' @export
saveDimReduction <- function(results, method, outdir=NULL, name=NULL,
    prefix=NULL){
    if ((is.null(name) || is.null(outdir)) && is.null(prefix)) {
        stop("Either name and outdir, or prefix have to be specified for saving
             dimensionality reduction results")
    }
    if (is.null(prefix)) {
        prefix <- paste(outdir, "/", name, "_", sep="")
    }
    saveRDS(results, paste(prefix, method, "_results.rds",
                                         sep=""))
    write.table(results$Yred, paste(prefix,  method,
                                              "_reducedY.csv", sep=""),
                quote=FALSE, row.names=TRUE, col.names=NA, sep=",")
    write.table(results$M, paste(prefix,  method, "_TandC.csv",
                         sep=""),
                quote=FALSE, row.names=TRUE, col.names=NA, sep=",")
}

#' Save results from analyseDimreduction
write_medians <- function(threshold) {
    sapply(threshold, function(thr) {
     write.table(paste(medians_maxcor$component[medians_maxcor$threshold >= 
                                                   thr],
                      collapse = ","),  
                paste(output ,"/ComponentsPassingThr", thr, "_", method, 
                      ".csv", sep=""), 
                col.names=FALSE, row.names=FALSE, quote=FALSE)
    })
}
