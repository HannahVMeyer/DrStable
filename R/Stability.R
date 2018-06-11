#' Estimate stability
#' 
#' @param dr [list]
#' @param threshold [double]
#' @param verbose [logical] If set, progress messages are printed to standard 
#' out
#' @return named list
#' @export
estimateStability <- function(dr, threshold, verbose=FALSE) {
    comparison  <- resultsComparison(dr, verbose=verbose)
    names(comparison) <- names(dr)
    formated <- formatComparison(comparison)
    maxcor <- formated$maxcor
    rownames(maxcor) <- unlist(sapply(1:(length(comparison)-1), function(i) {
        sapply(i:(length(comparison)-1), function(j, i) {
            paste(names(comparison)[i], names(comparison)[j+1],sep="_")
        }, i=i)
    }))
    maxcor_long <- reshape2::melt(maxcor)
    colnames(maxcor_long) <- c("comparison", "component", "value")
    medians_maxcor <- stats::aggregate(abs(value) ~ component, maxcor_long,
                                       median)
    colnames(medians_maxcor)[2] <- "abs_correlation" 
    medians_maxcor$threshold <- 
        sapply(medians_maxcor$abs_correlation, function(m) {
            # index of last TRUE
            threshold[length(which(sapply(threshold, function(thr) m >= thr )))]
        })
    maxcor_long$threshold <- sapply(maxcor_long$component, function(s) {
        medians_maxcor$threshold[medians_maxcor$component %in% s]
    })
    
    pass <- corrPass(formated$maxcor, threshold)
    colnames(pass) <- threshold
    rownames(pass) <- rownames(maxcor)
    pass_long <- reshape2::melt(pass)
    colnames(pass_long) <- c("comparison", "threshold","value")
    medians_pass <- stats::aggregate(value ~ threshold, pass_long, median)
    return(list(comparison=comparison, maxcor=maxcor, maxcor_long=maxcor_long, 
                medians_maxcor=medians_maxcor, pass=pass, pass_long=pass_long, 
                medians_pass=medians_pass, p_maxcor=p_maxcor, p_pass=p_pass ))
}

#' Apply compareSetup and format results
#' 
#' @param data_list [list]
#' @param verbose [logical] If set, progress messages are printed to standard 
#' out
#' 
#' @return named list with the size of the overlapping samples between samples 
#' sets (size), the reordered correlation coefficients (reorder), the maximum 
#' correlation (maxcor), the original column-wise correlation coefficients 
#' (componentxcomponent_r) and the reordered column-wise correlation
#' coefficients (componentxcomponent_r_reordered)
resultsComparison <- function(data_list, verbose=FALSE) {
    lapply(seq_along(data_list), function(perm1) {
        tmp <- lapply(seq_along(data_list), function(perm2, perm1) {
            if (perm1 <  perm2) {
                vmessage(c("Comparing list element", perm1, "with element", 
                           perm2), verbose=verbose)
                compareSets(data_list[[perm2]], data_list[[perm1]])
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
}

#' Compare dimensionality reduction across subsets
#' 
#' @param set1 [M1 x D] matrix with M1 samples and D dimensionality reduced
#' data P
#' @param set2 [M2 x D] matrix with M2 samples and D dimensionality reduced
#' data P
#' @param verbose [logical] If set, progress messages are printed to standard 
#' out
#' @return named list
compareSets <- function(set1, set2, verbose=FALSE) {
    common_subset <- findIntersect(set1, set2)
    size_subset <- dim(common_subset[[1]])[2]
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
#' @param comparison_list output [list] from resulsComparison containing with 
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
    maxcor <- do.call(rbind, sapply(sapply(maxcor,rmnulls), function(l) 
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
#' @param set1 [M1 x D] matrix with M1 samples and D dimensionality reduced
#' data P
#' @param set2 [M2 x D] matrix with M2 samples and D dimensionality reduced
#' data P
#' @return list with set1 and set2 filtered for overlapping samples and ordered
#' to the order of set1
findIntersect <- function(set1, set2) {
	set2 <- set2[,which(colnames(set2) %in% colnames(set1))]
	set1 <- set1[,which(colnames(set1) %in% colnames(set2))]
	set2 <- set2[,match(colnames(set1),colnames(set2))]
 	return(list(set1, set2))
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
	apply(set1, 1, function(fac1) {
		facxfac <- apply(set2, 1, function(fac2, fac1) {
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
	
	mat_in = matrix(0, nr=nrow(mat), nc=ncol(mat))
	iteration=1

	while(!sum(diag(mat)) == 0)  {
		highest_diag_index_memory=c()
		for (i in 1: nrow(mat)){
            vmessage(c("Iteration:", iteration, "i:", i, "\n"))
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
                exclude_from_max <- unique(c(index_max_diag, 
										which(diag_used == 1)))
			    max_col_index_max_diag <- which(abs(mat[index_max_diag,]) == 
											max(abs(mat[index_max_diag,
												-exclude_from_max])))
			    max_row_index_max_diag <- which(abs(mat[,index_max_diag]) == 
											max(abs(mat[-exclude_from_max,
												index_max_diag])))
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
				if ( abs(mat[ index_max_diag,max_col_index_max_diag]) >= 
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


