`simulate.gamm` <- function(object, nsim = 1, seed = NULL, newdata,
                            freq = FALSE, unconditional = FALSE, ...) {
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        runif(1)
    if (is.null(seed))
        RNGstate <- get(".Random.seed", envir = .GlobalEnv)
    else {
        R.seed <- get(".Random.seed", envir = .GlobalEnv)
        set.seed(seed)
        RNGstate <- structure(seed, kind = as.list(RNGkind()))
        on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
    }

    if (missing(newdata)) {
        newdata <- object$gam$model
    }

    ## random multivariate Gaussian
    ## From ?predict.gam in mgcv package
    ## Copyright Simon N Wood
    ## GPL >= 2
    rmvn <- function(n, mu, sig) { ## MVN random deviates
        L <- mroot(sig)
        m <- ncol(L)
        t(mu + L %*% matrix(rnorm(m * n), m, n))
    }

    Rbeta <- rmvn(n = nsim,
                  mu = coef(object$gam),
                  sig = vcov(object$gam, freq = freq,
                             unconditional = unconditional))
    Xp <- predict(object$gam, newdata = newdata, type = "lpmatrix")
    sims <- Xp %*% t(Rbeta)
    sims
}
