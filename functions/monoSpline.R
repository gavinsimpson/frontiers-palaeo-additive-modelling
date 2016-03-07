`monoSpline` <- function(depth, date, error = NULL, increasing = FALSE,
                         bs = "cr", k = 10, fx = FALSE, ...) {
    stopifnot(require("mgcv"))
    if(is.null(error)) {
        error <- rep(0, length(depth))
        W <- rep(1, length(depth))
    } else {
        W <- 1 / error
        W[1] <- W[2] * error[2]
    }
    df <- data.frame(Depth = depth, Date = date, error = error)
    mod <- gam(Date ~ s(Depth, k = k, bs = bs, fx = fx), data = df,
               fit = FALSE)
    sm <- smoothCon(s(Depth, k = k, bs = bs, fx = fx), data = df,
                    knots = NULL)[[1]]
    ## Fm are the constraints to enforce monotonicity
    Fm <- mono.con(sm$xp, up = increasing)
    G <- list(X = sm$X, C = matrix(0,0,0), sp = mod$sp, p = -sm$xp,
              y = df$Date, Ain = Fm$A, bin = Fm$b,
              S = sm$S, off = 0, w = W) ## w = df$Date*0+1
    ## fitted parameters
    p <- pcls(G)
    ## fitted values
    pm <- Predict.matrix(sm, data = data.frame(Depth = df$Depth))
    fit <- pm %*% p
    ## set up return object
    ## this needs all info for a predict() method to work, inc
    ## sm - smoothing matrix
    ## p  - a vector of fitted coefficients; the object returned from pcls
    ##
    ## also return fitted values
    out <- list(y = df$Date,
                fitted.values = drop(fit),
                coefficients = p,
                residuals = drop(df$Date - fit),
                sm = sm)
    class(out) <- "monoSpline"
    out
}

## predict S3 method
`predict.monoSpline` <- function(object, newdata, ...) {
    if(missing(newdata))
        yhat <- fitted(object)
    else {
        yhat <- Predict.matrix(object$sm, data = newdata) %*% coef(object)
    }
    yhat
}
