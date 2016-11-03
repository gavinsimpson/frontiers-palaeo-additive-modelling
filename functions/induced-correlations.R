
corFun <- function(V, t) {
    top <- V[t, ]
    bot <- sqrt(V[t,t] * diag(V))
    top / bot
}

inducedCor <- function(m, i) {
    ind <- rep(0,0)
    for (j in 1:length(m$smooth))
        ind <- c(ind, m$smooth[[j]]$first.para-1 + 1:m$smooth[[j]]$rank)
    mmat <- model.matrix(m)
    Z <- mmat[,ind]
    ZZt <- tcrossprod(Z)
    I <- diag(nrow(mmat))
    tmp <- capture.output(vcomps <- gam.vcomp(m))
    sigma <- vcomps[2,1]
    tau <- vcomps[1,1]
    V <- (sigma^2 * I) + (tau^2 * ZZt)
    corFun(V, i)
}

addInduced <- function(t, tp, mod) {
    cors <- unlist(lapply(tp, function(i, mod) inducedCor(mod, i), mod = mod),
                   use.names = FALSE)
    df <- data.frame(t  = rep(t, length(tp)),
                     tp = factor(rep(tp, each = length(t))),
                     correl = cors)
    df
}
