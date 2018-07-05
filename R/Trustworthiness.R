CalculateM <- function (k, original_rank, reduced_rank, N, random) {
    # define neighbourhoods based on ranks:
    # if xj part of neighbourhood -> 1 else 0
    # (set k+1 and exclude 1 due to rank(xi,xi)=1)
    if (! random) {
        original_neighbour <- t(apply(original_rank, 1, function(row)
            sapply(row, function(entry, k)
                if (entry <= k+1  && entry!=1) {return(1)}
                else {return(0)}, k=k)))
        reduced_neighbour <- t(apply(reduced_rank, 1, function(row)
            sapply(row, function(entry, k)
                if (entry <= k+1 && entry!=1) {return(1)}
                else {return(0)}, k=k)))
    } else {
        original_neighbour <- t(apply(original_rank, 1, function(row)
            sample(c(rep(1, k), rep(0, length(row)-k)))))
        reduced_neighbour <- t(apply(reduced_rank, 1, function(row)
            sample(c(rep(1, k), rep(0, length(row)-k)))))
    }
    
    # Data set U: find rank of samples xj which are within k in the reduced 
    # display but not in the original data set
    U <- sapply(seq(1, N, 1), function(x)
        original_rank[x,intersect(which(reduced_neighbour[x,] == 1),
                                  which(original_neighbour[x,] == 0))])
    
    # sum of all ranks(xi,xj) -k for all xj in U
    sum_ranksU <- sapply(U, function(element, k)
        if (length(element) == 0) {return(0)}
        else {return(sum(element) -k*length(element))}, k=k)
    
    # sum sum_ranks over all xi in i=1 to N
    sum_U <- sum(sum_ranksU)
    
    # Data set V: find samples xj which are within k in the original data
    # but not in the reduced data set
    V <- sapply(seq(1, N, 1), function(x)
        reduced_rank[x,intersect(which(original_neighbour[x,] == 1),
                                 which(reduced_neighbour[x,] == 0))])
    
    # sum of all ranks(xi,xj) -k for all xj in U
    sum_ranksV <- sapply(V, function(element, k)
        if (length(element) == 0) {return(0)}
        else {return(sum(element) -k*length(element))}, k=k)
    
    # sum sum_ranks over all xi in i=1 to N
    sum_V <- sum(sum_ranksV)
    
    # Calculate Trustworthiness M1
    M1 <- 1 - 2/(N*k*(2*N-3*k-1)) * sum_U
    
    # Calculate Continuity M2
    M2 <- 1 - 2/(N*k*(2*N-3*k-1)) * sum_V
    
    return(c(M1, M2))
}



TandC <- function(reduced, original, neighbours, random=FALSE){
    # compute distance matrices between samples (rows)
    # of high-dimensional data space (cols)
    if (class(original)=="dist") {
        original_dist <- as.matrix(original)
    } else {
        original_dist <- as.matrix(dist(original, method="euclidean"))
    }
    if (class(reduced)=="dist") {
        reduced_dist <- as.matrix(reduced)
    } else {
        reduced_dist <- as.matrix(dist(reduced, method="euclidean"))
    }
    
    N <- nrow(original_dist)
    
    # convert distance into ranks: rows xi, cols xj i.e. per row ranked
    # distances of xj with respective xi
    original_rank <- t(apply(original_dist, 1, function(x)
        rank(x, ties.method="first")
        ))
    reduced_rank <- t(apply(reduced_dist, 1, function(x)
        rank(x, ties.method="first")
    ))
    if (! random) {
        M <- sapply(neighbours, CalculateM, original_rank=original_rank,
                reduced_rank=reduced_rank, N=N, random=random)
    } else {
        M <- sapply(neighbours, CalculateM, original_rank=original_rank,
                    reduced_rank=original_rank, N=N, random=random)
        
    }
    rownames(M) <- c("T", "C")
    colnames(M) <- neighbours
    return(M)
}


