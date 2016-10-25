Mat<-
function (genos){
    X <- t(as.matrix(genos))
    n <- nrow(X)
    frac.missing <- apply(X, 2, function(x) {
        length(which(is.na(x)))/n
    })
    missing <- max(frac.missing) > 0
    freq <- apply(X + 1, 2, function(x) {
        mean(x, na.rm = missing)
    })/2
    MAF <- apply(rbind(freq, 1 - freq), 2, min)
    min.MAF <- 1/(2 * n)
    max.missing <- 1 - 1/(2 * n)
    markers <- which((MAF >= min.MAF) & (frac.missing <= max.missing))
    m <- length(markers)
    var.A <- 2 * mean(freq[markers] * (1 - freq[markers]))
    one <- matrix(1, n, 1)
    mono <- which(freq * (1 - freq) == 0)
    X[, mono] <- 2 * tcrossprod(one, matrix(freq[mono], length(mono), 1)) - 1
    freq.mat <- tcrossprod(one, matrix(freq[markers], m, 1))
    W <- X[, markers] + 1 - 2 * freq.mat
    if (missing){
        isna <- which(is.na(W))
        W[isna] <- 0
    }
    A <- tcrossprod(W)/var.A/m
    rownames(A) <- rownames(X)
    colnames(A) <- rownames(A)
    return(A)
}

