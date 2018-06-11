#' Save results of dimensionality reduction
#' 
#' @param results named list of results from dimensionality reduction; if 
#' multiple methods were specified, each named list entry containing:
#' Y_red: named list with dimensionality reduced phenotypes (reducedY) and 
#' object returned by specified dimensionality reduction method (results) with
#' additional output.
#' M: vector [double] with Trustworthiness and Continuity estimates for the
#' dimensionality reduction.
#' @param outdir path/to/directory [string] where results will be saved.
saveDimReduction <- function(results, outdir) {
    saveRDS(results$Y_red$results, paste(outdir, "/", m, "_results.rds", 
                                         sep=""))
    write.table(resutls$Y_red$reducedY, paste(outdir, ",",  m, "_reducedY.csv", 
                                    sep=""), 
                quote=FALSE, row.names=TRUE, col.names=NA, sep=",")
    write.table(results$M, paste(outdir, ",",  m, "_TandC.csv", 
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