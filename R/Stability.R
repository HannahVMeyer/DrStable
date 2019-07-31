#' Estimate stability
#' 
#' Takes the output of \code{\link{subsetDimReduction}} and finds
#' the stable low-dimension components across dimensionality reduced
#' subsets of the original data.
#' @param dr [list] either output of \code{\link{subsetDimReduction}} or
#' manually generated list with dimensionality reduction results on different
#' subsets of the same dataset. Each list entry is a matrix with row and column
#' names; row.names are crucial to allow for the comparison of the lower
#' dimensional embedding across common sample sets.
#' @param threshold [double] Threshold for stability between 0 and 1 or vector 
#' of thresholds between 0 and 1.
#' @param procrustes [logical] indication whether a procrustes transformation
#' to match the lower dimensional representations of common samples between
#' subsets should be performed.
#' @param verbose [logical] If set, progress messages are printed to standard
#' out
#' @return named list with [1] stability: list with one dataframe per threshold,
#' where each dataframe contains the stability estimate of the low-dimensional
#' components for that threshold value, [2] corr: data.frame with correlation of
#' the low-dimensional components for all pairwise comparisons and components.
#' @export
estimateStability <- function(dr, threshold, procrustes=FALSE, verbose=FALSE) {
    comparison  <- resultsComparison(dr, procrustes=procrustes, verbose=verbose)
    formated <- formatComparison(comparison)
    cor_matrix <- formated$maxcor
    rownames(cor_matrix) <- unlist(sapply(1:(length(comparison)-1), 
                                             function(i) {
        sapply(i:(length(comparison)-1), function(j, i) {
            paste(names(comparison)[i], names(comparison)[j+1],sep="_")
        }, i=i)
    }))
    
    cor_df <- reshape2::melt(cor_matrix)
    colnames(cor_df) <- c("comparison", "component", "correlation")
    cor_df$abs_correlation <- abs(cor_df$correlation)
    
    unique_components <- data.frame(component=unique(cor_df$component), 
                                    threshold=threshold)
    if (any(cor_df$abs_correlation > threshold)) {
        stability_count <- reshape2::acast(cor_df, component~., 
                                           value.var="abs_correlation",
                                           subset = plyr::.(cor_df$abs_correlation > 
                                                            threshold),
                                           length)
        stability_norm <- as.numeric(stability_count)/
            length(levels(cor_df$comparison))
        stability_comp <- 
            data.frame(stability=stability_norm,
                       component=as.numeric(rownames(stability_count)))
    } else {
        stability_comp <- data.frame(stability=rep(0, nrow(unique_components)),
                                     component=unique_components$component)
    }
    stability_all <- merge(stability_comp, unique_components, by="component",
                           all.y=TRUE)
    stability_all$stability[is.na(stability_all$stability)] <- 0
    return(list(stability=stability_all, corr=cor_df))
}

#' Compute median correlation of low-dimensional components across subsamples
#'
#' Takes the output of \code{\link{estimateStability}} and returns the median
#' correlation per component for the specified threshold.
#' @param es output [list] of \code{\link{estimateStability}} with named entries
#' [1] stability: list with one dataframe per threshold,
#' where each dataframe contains the stability estimate of the low-dimensional
#' components for that threshold value, [2] corr: data.frame with correlation of
#' the low-dimensional components for all pairwise comparisons and components.
#' @param threshold threshold [double] to mark median correlation as
#' above/below.
#' @return Dataframe with median correlation for each component.
medianCorr <- function(es, threshold) {
    medians_corr <- stats::aggregate(abs_correlation ~ component, es$corr,
                                   median)
    medians_corr$threshold <- 
        sapply(medians_corr$abs_correlation, function(m) {
            # index of last TRUE
            iT <- length(which(m >= threshold))
            if (iT == 0) return(paste("<", min(threshold),sep=""))
            if (iT != 0) return(threshold[iT])
        })
    return(medians_corr)
}

#' Apply compareSetup and format results
#'
#' @param data_list [list] output of \code{\link{subsetDimReduction}}. Each list
#' entry is a matrix with row and column names; row.names are crucial to allow
#' for the comparison of the lower dimensional embedding across common sample
#' sets.
#' @param procrustes [logical] indication whether a procrustes transformation
#' to match the lower dimensional representations of common samples between
#' subsets should be performed.
#' @param verbose [logical] If set, progress messages are printed to standard
#' out.
#' @return named list with the size of the overlapping samples between samples
#' sets (size), the reordered correlation coefficients (reorder), the maximum
#' correlation (maxcor), the original column-wise correlation coefficients
#' (componentxcomponent_r) and the reordered column-wise correlation
#' coefficients (componentxcomponent_r_reordered)
resultsComparison <- function(data_list, procrustes, verbose=FALSE) {
    if (is.null(names(data_list))) {
        names(data_list) <- paste("cv", 1:length(data_list), sep="")
    }
    tmp <- lapply(seq_along(data_list), function(perm1) {
        tmp <- lapply(seq_along(data_list), function(perm2, perm1) {
            if (perm1 <  perm2) {
                vmessage(c("Comparing list element", perm1, "with element",
                           perm2), verbose=verbose)
                compareSets(data_list[[perm2]]$Yred, data_list[[perm1]]$Yred,
                            procrustes=procrustes)
            }
        }, perm1=perm1)
        names(tmp) <- names(data_list)
        size <- lapply(tmp, function(perm) return(perm$size))
        reorder_correlation <- lapply(tmp, function(perm) 
            return(perm$reorder_correlation))
        reorder <- lapply(reorder_correlation, function(perm) 
            return(perm$reorder))
        maxcor <- lapply(reorder_correlation, function(perm) 
            return(perm$maxcor))
        componentxcomponent_r <- lapply(tmp, function(perm) 
            return(perm$componentxcomponent_r))
        componentxcomponent_r_reordered <- lapply(tmp, function(perm) 
            return(perm$componentxcomponent_r_reordered))
        return(list(size=size, reorder=reorder,
                    maxcor=maxcor, componentxcomponent_r=componentxcomponent_r, 
                    componentxcomponent_r_reordered=
                        componentxcomponent_r_reordered))
    })
    names(tmp) <- paste("cv", 1:length(data_list), sep="")
    return(tmp)
}

#' Compare dimensionality reduction across subsets
#'
#' @param set1 M1 x D [matrix] with M1 samples and D dimensionality reduced
#' data P.
#' @param set2 M2 x D [matrix] with M2 samples and D dimensionality reduced
#' data P.
#' @param procrustes [logical] indication whether a procrustes transformation
#' to match the lower dimensional representations of common samples between
#' set1 and set2 should be performed.
#' @param verbose [logical] If set, progress messages are printed to standard
#' out.
#' @return named list
compareSets <- function(set1, set2, verbose=FALSE, procrustes=FALSE) {
    common_subset <- findIntersect(set1, set2)
    size_subset <- dim(common_subset[[1]])[1]
    if (procrustes) {
        common_subset[[2]] <- MCMCpack::procrustes(common_subset[[1]],
                                                   common_subset[[2]],
                             dilation=TRUE, translation=TRUE)$X.new
    }
    componentxcomponent <- correlateSets(common_subset[[1]],
                                         common_subset[[2]])
    componentxcomponent_r <- sapply(componentxcomponent, function(fac) fac$r)
    componentxcomponent_p <- sapply(componentxcomponent, function(fac) fac$p)
    reorder_correlation <- analyseCorrelation(componentxcomponent_r,
                                              verbose=verbose)
    componentxcomponent_r_reordered <-
        componentxcomponent_r[reorder_correlation$reorder[1,],]
    return(list(size=size_subset, reorder_correlation=reorder_correlation,
                componentxcomponent_r=componentxcomponent_r,
                componentxcomponent_p=componentxcomponent_p,
                componentxcomponent_r_reordered=
                    componentxcomponent_r_reordered))
}

#' Format the comparison results
#'
#' Extract relevant data and remove null-entries
#' @param comparison_list output [list] from resultsComparison with
#' the size of the overlapping samples between samples sets (size), the
#' reordered correlation coefficients (reorder), the maximum correlation
#' (maxcor), the original column-wise correlation coefficients 
#' (componentxcomponent_r) and the reordered column-wise correlation
#' coefficients (componentxcomponent_r_reordered)
#'
#' @return list with same entries as input list, formated to only contain unique
#' comparisons: size of the overlapping samples between samples
#' sets (size), the reordered correlation coefficients (reorder), the maximum
#' correlation (maxcor), the original column-wise correlation coefficients
#' (componentxcomponent_r) and the reordered column-wise correlation
#' coefficients (componentxcomponent_r_reordered)
#'
formatComparison <- function(comparison_list) {
    sample_sizes <- sapply(comparison_list, function(perm) return(perm$size))
    reorder  <- lapply(comparison_list, function(perm)
        return(perm$reorder))[-length(comparison_list)]
    maxcor <- lapply(comparison_list, function(perm)
        return(perm$maxcor))[-length(comparison_list)]
    componentxcomponent_r  <- lapply(comparison_list, function(perm)
        return(perm$componentxcomponent_r))[-length(comparison_list)]
    componentxcomponent_r_reordered  <- lapply(comparison_list, function(perm)
        return(perm$componentxcomponent_r_reordered))[-length(comparison_list)]

    # rm all nulls
    maxcor <- do.call(rbind, sapply(sapply(maxcor, rmnulls), function(l)
        do.call(rbind,l)))
    reorder <- lapply(reorder, rmnulls)
    componentxcomponent_r <- lapply(componentxcomponent_r, rmnulls)
    componentxcomponent_r_reordered <- lapply(componentxcomponent_r_reordered,
                                              rmnulls)
    return(list(sample_sizes= sample_sizes, maxcor=maxcor, reorder=reorder,
                componentxcomponent_r=componentxcomponent_r,
                componentxcomponent_r_reordered= 
                    componentxcomponent_r_reordered))
}

#' Find intersecting samples between two sets and order according to first
#' 
#' @param set1 M1 x D [matrix] with M1 samples and D dimensionality reduced
#' data P.
#' @param set2 M2 x D [matrix] with M2 samples and D dimensionality reduced
#' data P.
#' @return [list] with set1 and set2 filtered for overlapping samples and
#' ordered to the order of set1.
findIntersect <- function(set1, set2) {
    set2 <- set2[which(rownames(set2) %in% rownames(set1)),]
    set1 <- set1[which(rownames(set1) %in% rownames(set2)),]
    set2 <- set2[match(rownames(set1),rownames(set2)),]
    return(list(set1, set2))
}

#' Aligng datasets with procrustes analyses
#'
#' @param set1 M x D [matrix] with M samples and D dimensionality reduced
#' data P; the matrix to be transformed.
#' @param set2 M x D [matrix] with M samples and D dimensionality reduced
#' data P (with the same sample order as set1); the target matrix.
#' @param dilation [logical] indicating whether set1 should be dilated.
#' @param translation [logical] indicating whether set1 should be translated.
#' @return [matrix] that is the Procrustes transformed version of set1.
alignSets <- function(set1, set2, dilation=TRUE, translation=TRUE) {
    res <- MCMCpack::procrustes(set1, set2,
                         dilation=dilation, translation=translation)
    return(res$X.new)
}

#' Column-wise correlation of two datasets
#'
#' @param set1 [M x D] matrix with M samples and D dimensionality reduced
#' data P
#' @param set2 [M x D] matrix with M samples and D dimensionality reduced
#' data P, with the same sample order as set1
#' @param type type [string] of correlation to use, one of spearman or pearson
#' @return named list with correlation coefficient (r) and p-value of
#' correlation (p) for each column-column correlation
correlateSets <- function(set1, set2, type="spearman") {
    apply(set1, 2, function(fac1) {
        facxfac <- apply(set2, 2, function(fac2, fac1) {
            tmp <- Hmisc::rcorr(cbind(fac2, fac1), type=type)
            return(list(r=tmp$r[1,2], p=tmp$P[1,2]))
        }, fac1=fac1)
    facxfac_r <- sapply(facxfac, function(fac) fac$r)
    facxfac_p <- sapply(facxfac, function(fac) fac$p)
    return(list(r=facxfac_r, p=facxfac_p))
    })
}

#' Order results based on highest correlation
#'
#' @param mat [D x D] matrix with correlation coefficients from column-wise
#' comparison of D lower dimensions
#' @param verbose [logical] If set, progress messages are printed to standard
#' out
#' @return named list with ordering of the correlations (reorder) and the
#' maximum correlations (maxcor)
analyseCorrelation <- function(mat, verbose=FALSE) {
    reorder <- rbind(rep(0, nrow(mat)), rep(0, nrow(mat)))
    maxcor <- rep(0, nrow(mat))
    diag_used <- rep(0, nrow(mat))

    mat_in = matrix(0, nrow=nrow(mat), ncol=ncol(mat))
    iteration=1
    medianTopCorr <- median(order(abs(diag(mat)), decreasing=TRUE)[1:3])
    if (medianTopCorr > ceiling(ncol(mat)/2)) {
        mat <- apply(mat, 1, rev)
        mat <- apply(mat, 1, rev)
        reversed <- TRUE
    } else {
        reversed <- FALSE
    }
    while(!sum(diag(mat)) == 0)  {
        highest_diag_index_memory=c()
        for (i in 1: nrow(mat)){
            vmessage(c("Iteration:", iteration, "i:", i, "\n"), verbose=FALSE)
            # Find the diagonal element with the highest correlation
            if (i == 1) {
                index_max_diag <- which.max(abs(diag(mat)))
            } else if ( i == nrow(mat)) {
                index_max_diag <- c(1:nrow(mat))[!c(1:nrow(mat)) %in%
                                        highest_diag_index_memory]
            } else {
                max_diag <- max(abs(diag(mat[-highest_diag_index_memory,
                                -highest_diag_index_memory])))
                index_max_diag <- which(abs(diag(mat)) == max_diag)
            }
            if (length(index_max_diag) > 1) {
                index_max_diag = index_max_diag[1]
            }

            # find maximum column and row element for max diagonal element
            if (iteration == 1) {
                max_col_index_max_diag <- which.max(abs(mat[index_max_diag,]))
                max_row_index_max_diag <- which.max(abs(mat[,index_max_diag]))
            } else {
                exclude_from_max <-
                    unique(c(index_max_diag, which(diag_used == 1)))
                max_col_index_max_diag <-
                    which(abs(mat[index_max_diag,]) ==
                        max(abs(mat[index_max_diag, -exclude_from_max])))
                max_row_index_max_diag <-
                    which(abs(mat[,index_max_diag]) ==
                        max(abs(mat[-exclude_from_max, index_max_diag])))
            }

            if (length(max_col_index_max_diag) > 1) {
                max_col_index_max_diag = max_col_index_max_diag[1]
            }
            if (length(max_row_index_max_diag) > 1) {
                max_row_index_max_diag = max_row_index_max_diag[1]
            }

            # check if the highest diagonal element is the highest element for
            # that row/column
            if (abs(mat[ index_max_diag, index_max_diag]) >=
                            abs(mat[ index_max_diag,max_col_index_max_diag]) &&
                abs(mat[ index_max_diag, index_max_diag]) >=
                            abs(mat[max_row_index_max_diag, index_max_diag]) &&
                diag_used[index_max_diag] != 1) {
                reorder[,index_max_diag] = t(rep(index_max_diag,2))
                maxcor[index_max_diag] = mat[ index_max_diag, index_max_diag]

                mat[index_max_diag,] <- 0
                mat[,index_max_diag] <- 0
            } else  {
                if (abs(mat[ index_max_diag,max_col_index_max_diag]) >=
                        abs(mat[max_row_index_max_diag, index_max_diag])) {
                    reorder[1,index_max_diag] <- max_col_index_max_diag
                    reorder[2,index_max_diag] <- index_max_diag
                    maxcor[index_max_diag] <- mat[index_max_diag, 
                                                    max_col_index_max_diag]

                    mat[-max_col_index_max_diag,max_col_index_max_diag] <- 0
                    mat[index_max_diag,] <- 0
                    diag_used[max_col_index_max_diag] <- 1
                } else {
                    reorder[1,index_max_diag] <- max_row_index_max_diag
                    reorder[2,index_max_diag] <- index_max_diag
                    maxcor[index_max_diag] <- mat[max_row_index_max_diag,
                                                index_max_diag]

                    mat[max_row_index_max_diag, -max_row_index_max_diag] <- 0
                    mat[,index_max_diag] <- 0
                    diag_used[max_row_index_max_diag] <- 1
                }
            }
            if (all(diag(mat) == 0)){
                if (reversed) {
                    reorder <- t(apply(reorder, 1, rev))
                    maxcor <- rev(maxcor)
                }
                return(list(reorder=reorder, maxcor=maxcor))
            }
            highest_diag_index_memory <- c(highest_diag_index_memory,
                                            index_max_diag)
        }
        iteration <- iteration + 1
    }
    return(list(reorder=reorder, maxcor=maxcor))
}

#' Filter correlations based on sequence of thresholds
#'
#' @param corrmat [M x D] correlation matrix
#' @param threshold stability threshold [double]
#' @return vector with number of stable components depending on threshold
corrPass <- function(corrmat, threshold) {
                sapply(threshold, function(thr, corrmat) {
                    apply(corrmat, 1, function(perm, thr)
                        length(which(abs(perm) >= thr)), thr=thr)
                 }, corrmat=corrmat)
}
