## Implements the power exponential and matern covariance functions

`powerExp` <- function(h, phi, k) {
    exp(- abs(h / phi)^k)
}

`matern` <- function(h, phi, k) {
    absh <- abs(h / phi)
    addn <- ifelse(k == 1.5, 0,
                   ifelse(k == 2.5, (1/3)*absh^2,
                          ((2/5)*absh^2) + ((1/15) * absh^3)))
    exp(- absh) * (1 + absh + addn)
}
