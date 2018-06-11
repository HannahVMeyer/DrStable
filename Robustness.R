#################
### Libraries ###
#################
libloc="/homes/hannah/R/x86_64-pc-linux-gnu-library/3.3"

library("data.table", lib.loc=libloc)
library("Hmisc", lib.loc=libloc)
library("corrplot", lib.loc=libloc)
library("ggplot2", lib.loc=libloc)
library("gplots", lib.loc=libloc)
library ("grid",lib.loc=libloc)
library("gridExtra",lib.loc=libloc)
library("reshape2",lib.loc=libloc)
library("wesanderson")
library ("R.methodsS3")
library ("R.oo")
library ("R.utils")
library("methods")

#################
### functions ###
#################

### Find intersecting samples between two sets and order according to first
findIntersect <- function(set1, set2) {
	set2 <- set2[,which(colnames(set2) %in% colnames(set1))]
	set1 <- set1[,which(colnames(set1) %in% colnames(set2))]
	set2 <- set2[,match(colnames(set1),colnames(set2))]
 	return(list(set1, set2))
}

### correlate two sets factor by factor and retrun r and pvalue
correlateSets <- function(set1, set2) {
	apply(set1, 1, function(fac1) {
		facxfac <- apply(set2, 1, function(fac2, fac1) {
			tmp <- rcorr(cbind(fac2, fac1), type="spearman")
			return(list(r=tmp$r[1,2], p=tmp$P[1,2]))
		}, fac1=fac1)
	facxfac_r <- sapply(facxfac, function(fac) fac$r)
	facxfac_p <- sapply(facxfac, function(fac) fac$p)
	return(list(r=facxfac_r, p=facxfac_p))
	})
}

### sum of squares
ss <- function(data) {
    return(sum(data^2))
}

### inverse normalise
invnorm = function(x) {
    qnorm((rank(x, na.last="keep", ties.method="random") - 1.5)/sum(!is.na(x)))
}

### order results based on highest correlation
analyseCorrelation <- function(mat, verbose=FALSE) {
	reorder = rbind(rep(0,nrow(mat)), rep(0,nrow(mat)))
	maxcor = rep(0,nrow(mat))
    diag_used = rep(0,nrow(mat))
	
	mat_in = matrix(0, nr=nrow(mat), nc=ncol(mat))
	iteration=1

	while( !sum(diag(mat)) == 0)  {
		highest_diag_index_memory=c()
		for (i in 1: nrow(mat)){
            if (verbose ) {
				cat("iteration:", iteration, "i:",i, "\n")
			}
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

### compare permuations
compareSetup <- function(set1, set2, verbose=FALSE) {
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

### remove nulls from nested list
rmnulls <- function(x) {
    x[!vapply(x, is.null, logical(1))]
}

### apply compareSetup and format results
resultsComparison <- function(data_list, verbose=FALSE) {
	lapply(seq_along(data_list), function(perm1) {
		tmp <- lapply(seq_along(data_list), function(perm2, perm1) {
            if (perm1 <  perm2) {
                if (verbose) {
 			        cat("Comparing list element", perm1, "with element", perm2)
		        }
                compareSetup(data_list[[perm2]], data_list[[perm1]])
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

### Format the comparison results: extract relevant data and remove null-entries
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
    reorder <- lapply(reorder,rmnulls)
    componentxcomponent_r <- lapply(componentxcomponent_r, rmnulls)
    componentxcomponent_r_reordered <- lapply(componentxcomponent_r_reordered,
																	rmnulls)
    return(list(sample_sizes= sample_sizes,maxcor=maxcor, reorder=reorder, 
				componentxcomponent_r=componentxcomponent_r,  
				componentxcomponent_r_reordered= 
											componentxcomponent_r_reordered))
}

### Filter correlations based on sequence of thresholds
corrPass <- function(corrmat, threshold) {
                sapply(threshold, function(thr, corrmat) {
                    apply(corrmat, 1, function(perm, thr) 
						length(which(abs(perm) >= thr)), thr=thr)
                 }, corrmat=corrmat)
}


analyseDimreduction <- function(dr, output, threshold, verbose=FALSE) {
    comparison  <- resultsComparison(dr, verbose=verbose)
    names(comparison) <- names(dr)
    formated <- formatComparison(comparison)
    maxcor <- formated$maxcor
    rownames(maxcor) <- unlist(sapply(1:(length(comparison)-1), function(i) {
        sapply(i:(length(comparison)-1), function(j, i) {
            paste(names(comparison)[i], names(comparison)[j+1],sep="_")
        }, i=i)
    }))
    maxcor_long <- melt(maxcor)
    colnames(maxcor_long) <- c("comparison", "component","value")
    medians_maxcor <- aggregate(abs(value) ~ component, maxcor_long, median)
    colnames(medians_maxcor)[2] <- "abs_correlation" 
    medians_maxcor$threshold <- 
		sapply(medians_maxcor$abs_correlation, function(m) {
        # index of last TRUE
        threshold[length(which(sapply(threshold, function(thr)  m >= thr )))]
    })
    maxcor_long$threshold <- sapply(maxcor_long$component, function(s) {
		medians_maxcor$threshold[medians_maxcor$component %in% s]
	})

    write_medians <- sapply(threshold, function(thr) {
		write.table(paste(medians_maxcor$component[medians_maxcor$threshold >= 
																		thr], 
						collapse = ","),  
					paste(output ,"/ComponentsPassingThr", thr, "_", method, 
                          ".csv", sep=""), 
					col.names=FALSE, row.names=FALSE, quote=FALSE)
	})
  
    
    pass <- corrPass(formated$maxcor, threshold)
    colnames(pass) <- threshold
    rownames(pass) <- rownames(maxcor)
    pass_long <- melt(pass)
    colnames(pass_long) <- c("comparison", "threshold","value")
    medians_pass <- aggregate(value ~ threshold, pass_long, median)
    
    pdf(file=paste(output, "/", strftime(Sys.time(), "%Y%m%d"),
		"_CorrelationComponentsBoxplots_", method, ".pdf", sep=""), 
        onefile=TRUE, height=12, width=16, paper = "a4r")
    p_pass <- ggplot(pass_long, aes(x=as.factor(threshold), y=value))
    p_pass <- p_pass +
        geom_boxplot(outlier.colour=NA) + 
        geom_jitter(width = 0.2, size=0.5) +
        labs(x = "Correlation threshold for filtering", 
			 y = "Number of components") +
        geom_text(data = medians_pass, aes(x = as.factor(threshold), 
										   y = value, label= value), 
				  nudge_x =-0.7, color="red") +
        theme_bw()
    
    p_maxcor <- ggplot(maxcor_long, aes(x=as.factor(component), y=abs(value), 
										fill=as.factor(threshold)))
    p_maxcor <- p_maxcor + 
        geom_boxplot(outlier.colour=NA) +
        scale_fill_manual(values=
			color[(16-length(levels(as.factor(maxcor_long$threshold)))):16]) +
        guides(fill=guide_legend(title="Threshold")) +
        labs(x = "Components", y = "Spearman correlation") +
        theme_bw() +
        theme(axis.text.x = element_text(angle=90))
    
    p_pass
    p_maxcor
    dev.off()
    return(list(comparison=comparison, maxcor=maxcor, maxcor_long =maxcor_long, 
				medians_maxcor=medians_maxcor, pass=pass, pass_long=pass_long, 
				medians_pass=medians_pass, p_maxcor=p_maxcor, p_pass=p_pass ))
}


############
### data ###
############

### read command line arguments
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


################
### analysis ###
################

if (length(dr) > 1) {
	cat("Determine robustness of ", method, " for ", N, " samples")
	color=c(wes_palette(20, name="Zissou", type='continuous'))
	threshold=seq(0.0, 0.95, 0.05)
	results <- analyseDimreduction(dr, output=output, threshold=threshold)
    saveRDS(results, paste(output, "/", method, "_robustness.rds", sep=""))
} else {
	cat("Not enough datasets (#", length(dr), ") to compute robustness")
}


