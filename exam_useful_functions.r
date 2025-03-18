library(rrBLUP)
make_kinship = function(genos){
    A.mat(t(genos[,-c(1:3)]))
}

qvalue <- function(p) {
    smooth.df = 3
    if (min(p) < 0 || max(p) > 1) {
        print("ERROR: p-values not in valid range.")
        return(0)
    }
    lambda = seq(0, 0.9, 0.05)
    m <- length(p)
    pi0 <- rep(0, length(lambda))
    for (i in 1:length(lambda)) {
        pi0[i] <- mean(p >= lambda[i])/(1 - lambda[i])
    }
    spi0 <- smooth.spline(lambda, pi0, df = smooth.df)
    pi0 <- predict(spi0, x = max(lambda))$y
    pi0 <- min(pi0, 1)
    if (pi0 <= 0) {
        print("ERROR: The estimated pi0 <= 0. Check that you have valid p-values.")
        return(0)
    }
    u <- order(p)
    qvalue.rank <- function(x) {
        idx <- sort.list(x)
        fc <- factor(x)
        nl <- length(levels(fc))
        bin <- as.integer(fc)
        tbl <- tabulate(bin)
        cs <- cumsum(tbl)
        tbl <- rep(cs, tbl)
        tbl[idx] <- tbl
        return(tbl)
    }
    v <- qvalue.rank(p)
    qvalue <- pi0 * m * p/v
    qvalue[u[m]] <- min(qvalue[u[m]], 1)
    for (i in (m - 1):1) {
        qvalue[u[i]] <- min(qvalue[u[i]], qvalue[u[i + 1]], 
                            1)
    }
    return(qvalue)
}
manhattan <- function(input, fdr.level = 0.05) {
    input <- input[order(input[, 2], input[, 3]), ]
    chroms <- unique(input[, 2])
    n.chrom <- length(chroms)
    chrom.start <- rep(0, n.chrom)
    chrom.mid <- rep(0, n.chrom)
    if (n.chrom > 1) {
        for (i in 1:(n.chrom - 1)) {
            chrom.start[i + 1] <- chrom.start[i] + max(input[which(input[, 
                                                                         2] == chroms[i]), 3]) + 1
        }
    }
    x.max <- chrom.start[n.chrom] + max(input[which(input[, 
                                                          2] == chroms[n.chrom]), 3])
    plot(0, 0, type = "n", xlim = c(0, x.max), ylim = c(0, 
                                                        max(input[, 4]) + 1), ylab = "-log(p)", xlab = "Chromosome", 
         xaxt = "n")
    for (i in seq(1, n.chrom, by = 2)) {
        ix <- which(input[, 2] == chroms[i])
        chrom.mid[i] <- median(chrom.start[i] + input[ix, 
                                                      3])
        points(chrom.start[i] + input[ix, 3], input[ix, 4], 
               col = "dark blue", pch = 16)
    }
    if (n.chrom > 1) {
        for (i in seq(2, n.chrom, by = 2)) {
            ix <- which(input[, 2] == chroms[i])
            chrom.mid[i] <- median(chrom.start[i] + input[ix, 
                                                          3])
            points(chrom.start[i] + input[ix, 3], input[ix, 
                                                        4], col = "cornflowerblue", pch = 16)
        }
    }
    q.ans <- qvalue(10^-input[, 4])
    temp <- cbind(q.ans, input[, 4])
    temp <- temp[order(temp[, 1]), ]
    if (temp[1, 1] < fdr.level) {
        temp2 <- tapply(temp[, 2], temp[, 1], mean)
        qvals <- as.numeric(rownames(temp2))
        x <- which.min(abs(qvals - fdr.level))
        first <- max(1, x - 2)
        last <- min(x + 2, length(qvals))
        if ((last - first) < 4) {
            last <- first + 3
        }
        splin <- smooth.spline(x = qvals[first:last], y = temp2[first:last], 
                               df = 3)
        lines(x = c(0, x.max), y = rep(predict(splin, x = fdr.level)$y, 
                                       2), lty = 2)
    }
    axis(side = 1, at = chrom.mid, labels = chroms)
}
qq <- function(scores) {
    remove <- which(scores == 0)
    if (length(remove) > 0) {
        x <- sort(scores[-remove], decreasing = TRUE)
    }
    else {
        x <- sort(scores, decreasing = TRUE)
    }
    n <- length(x)
    unif.p <- -log10(ppoints(n))
    plot(unif.p, x, pch = 16, xlab = "Expected -log(p)", 
         ylab = "Observed -log(p)")
    lines(c(0, max(unif.p)), c(0, max(unif.p)), lty = 2)
}