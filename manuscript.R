## ----knitr-defaults, cache = FALSE, echo = FALSE, message = FALSE, warning = FALSE----
knitr::opts_chunk$set(comment=NA, fig.align = "center", out.width = "0.8\\linewidth",
                      echo = FALSE, message = FALSE, warning = FALSE, cache = TRUE)
knitr::knit_hooks$set(crop.plot = knitr::hook_pdfcrop)

## ----load-packages, cache = FALSE----------------------------------------
library("mgcv")
library("scam")
library("ggplot2")
library("cowplot")
library("grid")                         # for unit.pmax(), unit.list()
library("schoenberg")
library("tidyr")

## Default ggplot theme
theme_set(theme_bw())

## ----common-components, cache = FALSE------------------------------------
## source Small Water data
small <- readRDS("./data/small-water/small-water-isotope-data.rds")

## load braya so data set
braya <- read.table("./data/braya-so/DAndrea.2011.Lake Braya So.txt", skip = 84)
names(braya) <- c("Depth", "DepthUpper", "DepthLower", "Year", "YearYoung", "YearOld", "UK37")
braya <- transform(braya, sampleInterval = YearYoung - YearOld)

## plot labels
d15n_label <- expression(delta^{15}*N)
braya_ylabel <- expression(italic(U)[37]^{italic(k)})

## ----data-figure, fig = TRUE, fig.width = 8, fig.height = 5, fig.cap = "Example time series; a) Small Water bulk organic matter $\\delta^{15}\\text{N}$ time series on a ${}^{210}\\text{Pb}$ time scale, and b) Braya-Sø \\uk{} time series on a calibrated ${}^{14}\\text{C}$ time scale. The observations \\uk{} time series have been joined by lines purely as a visual aid to highlight potential trends."----
## Generate a plot of the data - used as a base later
small_plt <- ggplot(small, aes(x = Year, y = d15N)) +
    geom_point() +
    labs(y = d15n_label, x = "Year CE")

## Generate a plot of the data
braya_plt <- ggplot(braya, aes(x = Year, y = UK37)) +
    geom_line(colour = "grey") +
    geom_point() +
    labs(y = braya_ylabel, x = "Year CE")

plot_grid(small_plt, braya_plt, ncol = 1, labels = "auto", align = "hv", axis = "lr")

## ----polynomial-example--------------------------------------------------
p <- c(1,3,5,10)
N <- 300
newd <- with(small, data.frame(Year = seq(min(Year), max(Year), length = N)))
polyFun <- function(i, data = data) {
    lm(d15N ~ poly(Year, degree = i), data = data)
}
mods <- lapply(p, polyFun, data = small)
pred <- vapply(mods, predict, numeric(N), newdata = newd)
colnames(pred) <- p
newd <- cbind(newd, pred)
polyDat <- gather(newd, Degree, Fitted, - Year)
polyDat <- transform(polyDat, Degree = ordered(Degree, levels = p))

## ----polynomial-example-plot, fig.width = 8, fig.height = 2.5, fig.cap = "Linear models with various orders of polynomial of the covariate Year fitted using ordinary least squares to the $\\delta^{15}\\text{N}$ time series from Small Water. The degree of polynomial is indicated, with the degree 1 line equal to a simple linear regression model.", dependson = c(-1)----
small_plt + geom_line(data = polyDat, mapping = aes(x = Year, y = Fitted, colour = Degree)) +
    scale_color_brewer(name = "Degree", palette = "PuOr") +
    theme(legend.position = "right")

## ----basis-function-example-plot, fig.width = 8, fig.height = 5, fig.cap = "Basis functions for the time covariate and the Small Water $\\delta^{15}\\text{N}$ time series. A rank (size) 7 cubic regression spline (CRS) basis expansion is show in a), with knots, indicated by tick marks on the x-axis, spread evenly through the rang of the data. b) shows the same CRS basis functions weighted by the estimated coefficients $\\beta_j$, plus the resulting GAM trend line (black line drawn through the data). The grey points in both panels are the observed $\\delta^{15}\\text{N}$ values. c) A rank 7 thin plate regression spline basis for the same data."----
## set up
k <- 7
df <- with(small, data.frame(Year = seq(min(Year), max(Year), length = 200)))
knots <- with(small, list(Year = seq(min(Year), max(Year), length = k)))
sm <- smoothCon(s(Year, k = k, bs = "cr"), data = df, knots = knots)[[1]]$X
colnames(sm) <- levs <- paste0("F", seq_len(k))
basis <- gather(cbind(sm, df), Fun, Value, -Year)
basis <- transform(basis, Fun = factor(Fun, levels = levs))

sm2 <- smoothCon(s(Year, k = k, bs = "cr"), data = small, knots = knots)[[1]]$X
beta <- coef(lm(d15N ~ sm2 - 1, data = small))
scbasis <- sweep(sm, 2L, beta, FUN = "*")
colnames(scbasis) <- levs <- paste0("F", seq_len(k))
fit <- cbind(df, fitted = rowSums(scbasis))
scbasis <- gather(cbind(scbasis, df), Fun, Value, -Year)
scbasis <- transform(scbasis, Fun = factor(Fun, levels = levs))

ylims <- range(basis$Value, scbasis$Value, small$d15N)

p1 <- ggplot(basis, aes(x = Year, y = Value, group = Fun, colour = Fun)) +
    geom_path() +
    scale_x_continuous(breaks = knots$Year, labels = NULL, minor_breaks = NULL) +
    scale_y_continuous(limits = ylims) +
    scale_colour_discrete(name = "Basis Function") +
    theme(legend.position = "none") +
    geom_point(data = small, mapping = aes(x = Year, y = d15N), inherit.aes = FALSE, size = 2, colour = "grey70") +
    labs(y = d15n_label, x = "Year CE (Knots)")

p2 <- ggplot(scbasis, aes(x = Year, y = Value, group = Fun, colour = Fun)) +
    geom_path() +
    scale_x_continuous(breaks = knots$Year, labels = NULL, minor_breaks = NULL) +
    scale_y_continuous(limits = ylims) +
    scale_colour_discrete(name = "Basis Function") +
    theme(legend.position = "none") +
    geom_point(data = small, mapping = aes(x = Year, y = d15N), inherit.aes = FALSE, size = 2, colour = "grey70") +
    geom_line(data = fit, mapping = aes(x = Year, y = fitted), inherit.aes = FALSE,
              size = 0.75, colour = "black") +
    labs(y = d15n_label, x = "Year CE (Knots)")

tp <- smoothCon(s(Year, k = k, bs = "tp"), data = df)[[1]]$X
colnames(tp) <- levs <- paste0("F", seq_len(k))
tpbasis <- gather(cbind(tp, df), Fun, Value, -Year)
tpbasis <- transform(tpbasis, Fun = factor(Fun, levels = levs))

p3 <- ggplot(tpbasis, aes(x = Year, y = Value, group = Fun, colour = Fun)) +
    geom_path() +
    scale_colour_discrete(name = "Basis Function") +
    theme(legend.position = "none") +
    labs(y = d15n_label, x = "Year CE")

pbasis <- plot_grid(p1, p2, p3, ncol = 2, align = "hv", labels = "auto")
pbasis

## ----penalty-example-----------------------------------------------------
K <- 40
lambda <- c(10000, 1, 0.01, 0.00001)
N <- 300
newd <- with(small, data.frame(Year = seq(min(Year), max(Year), length = N)))
fits <- lapply(lambda, function(lambda) gam(d15N ~ s(Year, k = K, sp = lambda), data = small))
pred <- vapply(fits, predict, numeric(N), newdata = newd)
op <- options(scipen = 100)
colnames(pred) <- lambda
newd <- cbind(newd, pred)
lambdaDat <- gather(newd, Lambda, Fitted, - Year)
lambdaDat <- transform(lambdaDat, Lambda = factor(paste("lambda ==", as.character(Lambda)), levels = paste("lambda ==", as.character(lambda))))
options(op)

## ----penalty-example-plot, fig.width = 8, fig.height = 5, fig.cap = "The effect of the smoothness parameter $\\lambda$ on the resulting wiggliness of the estimated spline. Large values of $\\lambda$ penalize wiggliness strongly, resulting in smooth trends (upper row), while smaller values allow increasingly wiggly trends. The aim of automatic smoothness selection is to find an optimal value of $\\lambda$ that balances the fit of the model with model complexity to avoid overfitting.", dependson = -1----
op <- options(scipen = 100)
small_plt + geom_line(data = lambdaDat, mapping = aes(x = Year, y = Fitted, group = Lambda),
                      size = .75, colour = "#e66101") +
    facet_wrap( ~ Lambda, ncol = 2, labeller = label_parsed)
options(op)

## ----fit-small-water-gamm------------------------------------------------
mod <- gamm(d15N ~ s(Year, k = 15), data = small,
            correlation = corCAR1(form = ~ Year), method = "REML")

## estimate of phi and confidence interval
smallPhi <- intervals(mod$lme, which = "var-cov")$corStruct

## summary object for use in document
smallSumm <- summary(mod$gam)

## ----fit-braya-so-car1-and-gcv-models------------------------------------
## fit the car(1) model --- needs optim as this is not a stable fit!
braya.car1 <- gamm(UK37 ~ s(Year, k = 10), data = braya, correlation = corCAR1(form = ~ Year),
                   method = "REML",
		   control = list(niterEM = 0, optimMethod = "BFGS", opt = "optim"))
braya.gcv <- gam(UK37 ~ s(Year, k = 30), data = braya)

## car(1) parameter
## optima-estimated version of the model gives phi = 0.200016
##           lower      est. upper
## Phi 2.58185e-16 0.2000162     1
brayaPhi <- intervals(braya.car1$lme)$corStruct

## ----plot-fitted-models, fig.width = 8, fig.height = 5, fig.cap = "GAM-based trends fitted to the Small Water $\\delta^{15}\\text{N}$ (a) and Braya-Sø \\uk{} (b) time series. The shaded bands surrounding the estimated trends are approximate 95\\% across-the-function confidence intervals. For the \\uk{} series, two models are shown; the orange fit is the result of a GAM with a continuous-time AR(1) process estimated using REML smoothness selection, while the blue fit is that of a simple GAM with GCV-based smoothness selection. The REML-based fit significantly oversmooths the \\uk{} time series."----
N <- 300                                # number of points at which to evaluate the splines
## Predict from the fitted model
newYear <- with(small, data.frame(Year = seq(min(Year), max(Year), length.out = 200)))
newYear <- cbind(newYear, data.frame(predict(mod$gam, newYear, se.fit = TRUE)))
newYear <- transform(newYear, upper = fit + (2 * se.fit), lower = fit - (2 * se.fit))

## Plot simulated trends
small_fitted <- ggplot(newYear, aes(x = Year, y = fit)) +
    geom_ribbon(aes(ymin = lower, ymax = upper, x = Year), alpha = 0.2, inherit.aes = FALSE, fill = "black") +
    geom_point(data = small, mapping = aes(x = Year, y = d15N), inherit.aes = FALSE) +
    geom_line() +
    labs(y = d15n_label, x = "Year CE")

## ggplot with data and fitted spline, then resids vs time in second panel
newBraya <- with(braya, data.frame(Year = seq(min(Year), max(Year), length.out = N)))
newBraya <- cbind(newBraya, data.frame(predict(braya.car1$gam, newBraya, se.fit = TRUE)))
crit.t <- qt(0.975, df = df.residual(braya.car1$gam))
newBraya <- transform(newBraya,
                      upper = fit + (crit.t * se.fit),
                      lower = fit - (crit.t * se.fit))
## add GAM GCV results
fit.gcv <- predict(braya.gcv, newdata = newBraya, se.fit = TRUE)
newBraya <- rbind(newBraya, newBraya)      # extend newBraya to take GCV results
newBraya[seq(N+1, length.out = N, by = 1), ]$fit <- fit.gcv$fit
newBraya[seq(N+1, length.out = N, by = 1), ]$upper <-
    fit.gcv$fit + (qt(0.975, df.residual(braya.gcv)) * fit.gcv$se.fit)
newBraya[seq(N+1, length.out = N, by = 1), ]$lower <-
    fit.gcv$fit - (qt(0.975, df.residual(braya.gcv)) * fit.gcv$se.fit)
newBraya <- transform(newBraya, Method = rep(c("GAMM (CAR(1))", "GAM (GCV)"), each = N))

## plot CAR(1) and GCV fits
braya_fitted <- ggplot(braya, aes(y = UK37, x = Year)) +
    geom_point() +
    geom_ribbon(data = newBraya,
                mapping = aes(x = Year, ymax = upper, ymin = lower, fill = Method),
                alpha = 0.3, inherit.aes = FALSE) +
    geom_line(data = newBraya, mapping = aes(y = fit, x = Year, colour = Method)) +
    labs(y = braya_ylabel, x = "Year CE") +
    scale_color_manual(values = c("#5e3c99", "#e66101")) + ## values = viridis(11)[c(2, 6)]) +
    scale_fill_manual(values = c("#5e3c99", "#e66101")) + ## values = viridis(11)[c(2, 6)]) +
    theme(legend.position = "right")

plot_grid(small_fitted, braya_fitted, ncol = 1, labels = "auto", align = "hv", axis = "lr")

## ----trace-smoothness-parameters, fig.width = 8, fig.height = 2.5, fig.cap = "GCV and REML scores as a function of the smoothness parameter $\\lambda$. From left to right, GAMs were estimated using GCV and REML smoothness selection, and REML using a basis dimension of 40 and observational weights to account for heterogeneity in the \\uk{} times series. The selected value of $\\lambda$ for each model is indicated by the vertical grey line."----
## Generate GCV and REML traces
lambda <- 1e-8
gcv <- reml <- reml2 <- numeric(length = 100)

for (i in seq_along(gcv)) {
    mGCV <- gam(UK37 ~ s(Year, k = 30), data = braya, method = "GCV.Cp", sp = lambda)
    mREML <- gam(UK37 ~ s(Year, k = 30), data = braya, method = "REML", sp = lambda)
    mREML2 <- gam(UK37 ~ s(Year, k = 40), data = braya, method = "REML", sp = lambda,
                  weights = sampleInterval / mean(sampleInterval))
    gcv[i] <- mGCV$gcv.ubre
    reml[i] <- mREML$gcv.ubre
    reml2[i] <- mREML2$gcv.ubre
    lambda <- lambda * 1.5
}

## GCV and REML score traces
result <- data.frame(lambda = 1e-8 * 1.5^{0:99}, score = c(gcv, reml, reml2),
                     criterion = rep(c("GCV", "REML", "REML (k=40, weights)"), each = 100))

## Full fits to get estimated models
fullGCV  <- gam(UK37 ~ s(Year, k = 30), data = braya, method = "GCV.Cp")
fullREML <- gam(UK37 ~ s(Year, k = 30), data = braya, method = "REML")
fullREML2 <- gam(UK37 ~ s(Year, k = 40), data = braya, method = "REML",
                 weights = sampleInterval / mean(sampleInterval))
obsSP <- data.frame(lambda = c(fullGCV$sp, fullREML$sp, fullREML2$sp), criterion = c("GCV", "REML", "REML (k=40, weights)"))

## summary object for use in document
brayaSumm <- summary(fullREML2)

## Plot
trace.plt <- ggplot(result, aes(x = lambda, y = score)) +
    geom_line() +
    facet_wrap(~ criterion, scales = "free_y") + scale_x_log10() +
    xlab(expression(lambda)) +
    ylab("Score") +
    geom_vline(aes(xintercept = lambda), obsSP, colour = "grey")
trace.plt

## ----posterior-simulation, fig.height = 5, fig.width = 8, fig.cap = "Estimated trends (thick black lines) and 20 random draws (grey lines) from the posterior distribution of the GAM fitted to the Small Water $\\delta^{15}\\text{N}$ (a) and Braya-Sø \\uk{} (b) time series."----
set.seed(1)
nsim <- 20
take <- 20
sims <- simulate(mod, nsim = nsim, newdata = newYear, unconditional = TRUE)
randSims <- sims[, sample(nsim, take)]
colnames(randSims) <- paste0("sim", seq_len(take))
randSims <- setNames(stack(as.data.frame(randSims)), c("simulated", "run"))
randSims <- transform(randSims, Year = rep(newYear$Year, take),
                      simulated = simulated)

## Plot simulated trends
smallSim.plt <- ggplot(newYear, aes(x = Year, y = fit)) +
    geom_line(data = randSims, mapping = aes(y = simulated, x = Year, group = run),
              colour = "grey80") +
    geom_line(lwd = 1) +
    labs(y = d15n_label, x = "Year CE")

## posterior simulation
## need to reset-up newBraya
newBraya <- with(braya, data.frame(Year = seq(min(Year), max(Year), length.out = N)))
brayaREML2 <- cbind(newBraya, data.frame(predict(fullREML2, newBraya, se.fit = TRUE)))
## simulate
set.seed(1)
sims2 <- simulate(fullREML2, nsim = nsim, newdata = newBraya, unconditional = TRUE)
randSims2 <- sims2[, sample(nsim, take)]
colnames(randSims2) <- paste0("sim", seq_len(take))
randSims2 <- setNames(stack(as.data.frame(randSims2)), c("simulated", "run"))
randSims2 <- transform(randSims2, Year = rep(newBraya$Year, take),
                       simulated = simulated)

brayaSim.plt <- ggplot(brayaREML2, aes(x = Year, y = fit)) +
    geom_line(data = randSims2, mapping = aes(y = simulated, x = Year, group = run),
              colour = "grey80") +
    geom_line(lwd = 1) +
    labs(y = braya_ylabel, x = "Year CE")

plot_grid(smallSim.plt, brayaSim.plt, ncol = 1, labels = "auto", align = "hv", axis = "lr")

## ----compare-intervals, fig.height = 5, fig.width = 8, fig.cap = "95\\% simultaneous confidence intervals (light grey bands) and across-the-function confidence intervals (dark grey bands) on the estimated trends (black lines) for the Small Water $\\delta^{15}\\text{N}$ (a) and Braya-Sø \\uk{} (b) time series."----
## small water
sw.cint <- confint(mod, parm = "Year", newdata = newYear, type = "confidence")
sw.sint <- confint(mod, parm = "Year", newdata = newYear, type = "simultaneous")

## braya so
bs.cint <- confint(fullREML2, parm = "Year", newdata = newBraya, type = "confidence")
bs.sint <- confint(fullREML2, parm = "Year", newdata = newBraya, type = "simultaneous")

smallInt.plt <- ggplot(sw.cint, aes(x = Year, y = est)) +
    geom_ribbon(data = sw.sint, mapping = aes(ymin = lower, ymax = upper, x = Year),
                fill = "grey80", inherit.aes = FALSE) +
    geom_ribbon(mapping = aes(ymin = lower, ymax = upper, x = Year),
                fill = "grey60", inherit.aes = FALSE) +
    geom_line(lwd = 1) +
    labs(y = d15n_label, x = "Year CE")

brayaInt.plt <- ggplot(bs.cint, aes(x = Year, y = est)) +
    geom_ribbon(data = bs.sint, mapping = aes(ymin = lower, ymax = upper, x = Year),
                fill = "grey80", inherit.aes = FALSE) +
    geom_ribbon(mapping = aes(ymin = lower, ymax = upper, x = Year),
                fill = "grey60", inherit.aes = FALSE) +
    geom_line(lwd = 1) +
    labs(y = braya_ylabel, x = "Year CE")

plot_grid(smallInt.plt, brayaInt.plt, ncol = 1, labels = "auto", align = "hv", axis = "lr")

## ----derivatives, fig.height = 5, fig.width = 8, fig.cap = "Estimated first derivatives (black lines) and 95\\% simultaneous confidence intervals of the GAM trends fitted to the Small Water $\\delta^{15}\\text{N}$ (a) and Braya-Sø \\uk{} (b) time series. Where the simultaneous interval does not include 0, the models detect significant temporal change in the response."----
small.d <- fderiv(mod, newdata = newYear, n = N)
small.sint <- with(newYear, cbind(confint(small.d, nsim = nsim, type = "simultaneous"), Year = Year))

small_deriv_plt <- ggplot(small.sint, aes(x = Year, y = est)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "black") +
    geom_line() +
    labs(x = "Year CE", y = "First derivative")

braya.d <- fderiv(fullREML2, newdata = newBraya, n = N)
braya.sint <- with(newBraya, cbind(confint(braya.d, nsim = nsim, type = "simultaneous"), Year = Year))

braya_deriv_plt <- ggplot(braya.sint, aes(x = Year, y = est)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "black") +
    geom_line() +
    labs(x = "Year CE", y = "First derivative")

plot_grid(small_deriv_plt, braya_deriv_plt, ncol = 1, labels = "auto", align = "hv", axis = "lr")

## ----car1-plot, fig.height = 2.5, fig.width = 4, out.width = "0.5\\linewidth", fig.cap = "Estimated CAR(1) process from the GAM fitted to the Small Water $\\delta^{15}\\text{N}$ time series. $h(\\Delta_t, \\phi)$ is the correlation between residuals separated by $\\Delta_t$ years, where $\\hat{\\phi} = \\text{0.6}$. The shaded band is a 95\\% pointwise confidence interval on the estimated correlation $h$."----
## plot CAR(1) process
maxS <- with(small, diff(range(Year))) ## too large, truncate to 50
S <- seq(0, 50, length = 100)

car1 <- setNames(as.data.frame(t(outer(smallPhi, S, FUN = `^`)[1, , ])),
                 c("Lower","Correlation","Upper"))
car1 <- transform(car1, S = S)

car1Plt <- ggplot(car1, aes(x = S, y = Correlation)) +
    geom_ribbon(aes(ymax = Upper, ymin = Lower),
                fill = "black", alpha = 0.2) +
    geom_line() +
    ylab(expression(italic(h) * (list(Delta[t], varphi)))) +
    xlab(expression(Delta[t] ~ (years)))
car1Plt

## ----gp-correlation-functions, fig.width = 8, fig.height = 5, fig.cap = "Power exponential (a) and Matérn (b) correlation functions for observation separation distance $h$. Two values of the effective range parameter ($\\phi$)) are shown for each function. For the power exponential function, $\\kappa$ is the power in the power exponential function. For the Matérn correlation function, $\\kappa$ distinguishes the member of the Matérn family."----
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

pExpdf <- expand.grid(h = seq(0, 7, length = 100), phi = c(1,2), k = c(0.5, 1, 2))
pExpdf <- transform(pExpdf, correlation = powerExp(h, phi, k), kfac = factor(k))

pexp.plt <- ggplot(pExpdf, aes(x = h, y = correlation, group = kfac, colour = kfac)) +
    geom_line(size = 1) +
    facet_wrap( ~ phi, labeller = label_bquote(phi == .(phi))) +
    ylab(expression(rho)) +
    xlab(expression(italic(h))) +
    ## scale_color_viridis(option = "viridis", discrete = TRUE) +
    ## scale_colour_manual(values = viridis(11)[c(2,5,10)]) +
    scale_colour_manual(name = expression(kappa), values = c("#e66101","#fdb863","#5e3c99"))

pMaternDf <- expand.grid(h = seq(0, 10, length = 100),
                         phi = c(1,2),
                         k = c(1.5, 2.5, 3.5))
pMaternDf <- transform(pMaternDf, correlation = matern(h, phi, k),
                       kfac = factor(k))

pmat.plt <- ggplot(pMaternDf, aes(x = h, y = correlation, group = kfac, colour = kfac)) +
    geom_line(size = 1) +
    facet_wrap( ~ phi, labeller = label_bquote(phi == .(phi))) +
    ylab(expression(rho)) +
    xlab(expression(italic(h))) +
    scale_colour_manual(name = expression(kappa), values = c("#e66101","#fdb863","#5e3c99"))

correlFuns <- plot_grid(pexp.plt, pmat.plt, ncol = 1, align = "hv",
                        labels = c("a","b"), axis = "lr")
correlFuns

## ----profile-range-parameter-in-gp-smooth-model--------------------------
nn <- 200                               # number of points at which to evaluate profile likelihood
dseq <- seq(15, 500, length.out = nn)  # effective ranges to fit at
## Mat <- SExp <- vector("list", length = nn)     # object to hold model fits
Mat <- SEx <- numeric(length = nn)     # object to hold model fits
for (i in seq_along(dseq)) {            # iterate over dseq, fit GP GAM w Matérn covariance
    Mat[i] <- gam(UK37 ~ s(Year, k = 45, bs = "gp", m = c(3, dseq[i])),
                  weights = sampleInterval / mean(sampleInterval),
                  data = braya, method = "REML", family = gaussian())[["gcv.ubre"]]
    SEx[i] <- gam(UK37 ~ s(Year, k = 45, bs = "gp", m = c(2, dseq[i], 1)),
                  weights = sampleInterval / mean(sampleInterval),
                  data = braya, method = "REML", family = gaussian())[["gcv.ubre"]]
}
## extract the REML score into ggplot-friendly object
## reml.scr <- data.frame(effrange = dseq, reml = sapply(mods, `[[`, "gcv.ubre"), aic = sapply(mods, AIC))
reml.scr <- data.frame(cor = rep(c("Matérn","Exponential"), each = nn),
                       effrange = rep(dseq, 2),
                       reml = c(Mat, SEx))

## ----fit-example-gp-smooths, dependsons = "profile-range-parameter-in-gp-smooth-model"----
effRange1 <- 250      # sets effective range in years for Matérn correl
effRange2 <- with(subset(reml.scr, cor == "Matérn"), dseq[which.min(reml)])
effRange3 <- with(subset(reml.scr, cor == "Exponential"), dseq[which.min(reml)])
## Add time interval of each sample as prior weights; effectively averaging more time per sample
## gp1 <- gam(UK37 ~ s(Year, k = 45, bs = "gp", m = c(3, effRange1)), data = braya,
##            method = "REML", weights = sampleInterval / mean(sampleInterval))
gp2 <- gam(UK37 ~ s(Year, k = 45, bs = "gp", m = c(3, effRange2)), data = braya,
           method = "REML", weights = sampleInterval / mean(sampleInterval))
gp3 <- gam(UK37 ~ s(Year, k = 45, bs = "gp", m = c(2, effRange3, 1)), data = braya,
           method = "REML", weights = sampleInterval / mean(sampleInterval))
newd <- with(braya, data.frame(Year = seq(min(Year), max(Year), length.out = 1000)))
## p.gp1 <- transform(newd, fitted = predict(gp1, newdata = newd, type = "response"),
##                    effRange = effRange1)
p.gp2 <- transform(newd, fitted = predict(gp2, newdata = newd, type = "response"),
                   effRange = round(effRange2))
p.gp3 <- transform(newd, fitted = predict(gp3, newdata = newd, type = "response"),
                   effRange = round(effRange3))
## pred <- rbind(p.gp1, p.gp2, p.gp3)
pred <- rbind(p.gp2, p.gp3)
pred <- transform(pred, effRange = factor(effRange),
                  cor = rep(c("Matérn", "Exponential"), each = nrow(newd)))

## ----gp-gam-detail-plot, fig.height = 5, fig.width = 8, fig.cap = "Gaussian process smooths fitted to the \\uk{} time series. REML score traces for GAMs fitted using power exponential ($\\kappa = 1$) or Matérn ($\\kappa = 1.5$) correlation functions as a function of the effective range parameter ($\\phi$) are shown (a). The optimal model for each function is that with the lowest REML score. b) shows the resulting trends estimated using the respective correlation function with the value of $\\phi$ set to the optimal value.", dependson="fit-example-gp-smooths"----
## profile-likelihood plot
proflik.plt <- ggplot(reml.scr, aes(x = effrange, y = reml, colour = cor)) +
    geom_line() +
    scale_colour_manual(name = "", values = c("#e66101","#5e3c99")) +
    labs(y = "REML score", x = expression(Effective ~ range ~ (varphi)))
## plot at two values of h
gp.plt2 <- ggplot(pred, aes(x = Year, y = fitted, colour = cor)) +
    geom_line() + theme(legend.position = "right") +
    geom_point(aes(x = Year, y = UK37), data = braya, inherit.aes = FALSE) +
    scale_colour_manual(name = "", values = c("#e66101","#5e3c99")) +
    labs(y = braya_ylabel, x = "Year CE")
plot_grid(proflik.plt, gp.plt2, ncol = 1, labels = c("a","b"), align = "hv", axis = "lr")

## ----braya-so-model-comparisons, fig.width = 8, fig.height = 2.5, fig.cap = "Comparison of trends estimated using i) adaptive smooth, ii) Gaussian process, and iii) thin plate regression spline bases for the \\uk{} time series.", dependson = "fit-example-gp-smooths"----
## model it using gam()
effRange <- effRange2

## Gaussian process, Matern, kappa = 1.5, weights as sampleInterval
mod_gp <- gam(UK37 ~ s(Year, k = 45, bs = "gp", m = c(3, effRange)), data = braya,
              method = "REML", weights = sampleInterval / mean(sampleInterval))

## Adaptive spline, weights as sampleInterval
mod_ad <- gam(UK37 ~ s(Year, k = 45, bs = "ad"), data = braya,
              method = "REML", weights = sampleInterval / mean(sampleInterval))

## TPRS, weights as sampleInterval, k = needs to be higher
mod_tprs <- gam(UK37 ~ s(Year, k = 45, bs = "tp"), data = braya, method = "REML",
                weights = sampleInterval / mean(sampleInterval))

## wrap this in a function that will return all the plots and derived objects
processGAM <- function(mod) {
    ## Predict from model
    N <- 500
    newYear <- with(braya, data.frame(Year = seq(min(Year), max(Year), length.out = N)))
    newYear <- cbind(newYear, data.frame(predict(mod, newYear, se.fit = TRUE)))
    
    out <- list(#plots = list(combinedPlt, plt.fit, plt.fit2, plt.sim, derivPlt),
                objects = newYear)
    out
}

plts_gp   <- processGAM(mod = mod_gp)      # Gaussian process smooth with weights
plts_ad   <- processGAM(mod = mod_ad)      # Adaptive smooth with weights
plts_tprs <- processGAM(mod = mod_tprs)    # TPRS with weights

pltData <- do.call("rbind", lapply(list(plts_gp, plts_ad, plts_tprs),
                                   `[[`, "objects"))
pltData <- transform(pltData, Model = rep(c("GP", "Adaptive", "TPRS"),
                              each = nrow(plts_gp$objects)))

allFits <- ggplot(pltData, aes(x = Year, y = fit)) +
    geom_point(aes(x = Year, y = UK37), data = braya) +
    geom_line(aes(colour = Model)) + labs(y = braya_ylabel, x = "Year") +
    theme(legend.position = "right") +
    scale_colour_manual(name = "", values = c("#e66101", "#fdb863", "#5e3c99"))
allFits

## ----small-scam-fit------------------------------------------------------
knots <- with(small, list(Year = seq(min(Year), max(Year), length = 14)))
mod <- gamm(d15N ~ s(Year, k = 15), data = small, method = "REML",
            correlation = corCAR1(form = ~ Year),
            knots = knots)

swAge <- read.csv("./data/small-water/small1-dating.csv")

## monotonic spline age-depth model
swAge$Error[1] <- 1.1
swAgeMod <- scam(Date ~ s(Depth, k = 5, bs = "mpd"), data = swAge, weights = 1 / swAge$Error, gamma = 1.4)

## predict from the age model for a smooth set of points in `Depth`
newAge <- with(swAge, data.frame(Depth = seq(min(Depth), max(Depth), length.out = 200)))
newAge <- transform(newAge, fitted = predict(swAgeMod, newdata = newAge, type = "response"))
newSims <- as.data.frame(simulate(swAgeMod, nsim = 25, newdata = newAge))
newSims <- cbind(Depth = newAge$Depth, newSims)
newSims <- gather(newSims, Simulation, Age, -Depth)

## simulate from age model; each column is a simulation
ageSims <- simulate(swAgeMod, nsim = 100, newdata = small, seed = 42)
ageSims <- as.data.frame(ageSims)

fitSWModels <- function(x, y, knots) {
    dat <- data.frame(d15N = y, Year = x)
    m <- gamm(d15N ~ s(Year, k = 15), data = dat, method = "REML",
              correlation = corCAR1(form = ~ Year), knots = knots)
}

predSWModels <- function(mod, newdata) {
    predict(mod$gam, newdata = newdata, type = "response")
}

simulateSWModels <- function(mod, newdata, nsim, seed = 42) {
    sims <- simulate(mod, nsim = nsim, newdata = newdata, seed = seed)
    as.vector(sims)
}

simTrendMods <- lapply(ageSims, fitSWModels, y = small$d15N, knots = knots)
simTrends <- lapply(simTrendMods, predSWModels, newdata = newYear)
simTrends <- data.frame(Year  = with(newYear, rep(Year, length(simTrends))),
                        Trend = unlist(simTrends),
                        Group = rep(seq_along(simTrends), times = lengths(simTrends)))

NSIM <- 50
simSimulate <- lapply(simTrendMods, simulateSWModels, newdata = newYear, nsim = NSIM, seed = 42)
simSimulate <- data.frame(Year  = with(newYear, rep(Year, times = NSIM * length(simSimulate))),
                          Trend = unlist(simSimulate),
                          Group = rep(seq_len(NSIM * length(simSimulate)), each = nrow(newYear)))

## ----small-scam-fit-plots, fig.width = 8, fig.height = 7.5, fig.cap = "Accounting for uncertainty in age estimates whilst fitting a smooth trend to the Small Water $\\delta^{15}\\text{N}$ time series. (a) Estimated age model using a monotonically-constrained spline fitted to ${}^{210}\\text{Pb}$ inferred ages for selected depths in the sediment core (red points). The uncertainty in the ${}^{210}\\text{Pb}$ inferred age is show by the red vertical bars. The fitted age model is illustrated by the solid black line. The faint grey lines are 25 random draws from the posterior distribution of the monotonically constrained GAM. The effect of age uncertainty on trend estimation is shown in b); for 100 simulations from the posterior distribution of the age model in a) a trend was estimated using a GAM with a thin plate regression spline basis and a CAR(1) process in the residuals. These trends are shown as grey lines. The combined effect of age model and fitted GAM uncertainty on the trends for the $\\delta^{15}\\text{N}$ time series is shown in c). The grey lines in c) are based on 50 random draws from the model posterior distribution for each of the 100 trends shown in b). For both b) and c) the black line shows the trend estimated assuming the ages of each sediment sample are known and fixed.", dependson = -1----
plt1 <- ggplot(swAge, aes(y = Date, x = Depth)) +
    geom_line(data = newSims, mapping = aes(y = Age, x = Depth, group = Simulation),
              alpha = 1, colour = "grey80") +
    geom_line(data = newAge, mapping = aes(y = fitted, x = Depth)) +
    geom_point(size = 1.5, colour = "red") +
    geom_errorbar(aes(ymin = Date - Error, ymax = Date + Error, width = 0), colour = "red") +
    labs(y = "Year CE", x = "Depth (cm)")

plt2 <- ggplot(simTrends, aes(x = Year, y = Trend, group = Group)) +
    geom_line(alpha = 0.1, colour = "grey80") +
    geom_line(data = newYear, mapping = aes(x = Year, y = fit), inherit.aes = FALSE) +
    geom_point(data = small, mapping = aes(x = Year, y = d15N), inherit.aes = FALSE, size = 0.7) +
    labs(x = "Year", y = d15n_label)

plt3 <- ggplot(simSimulate, aes(x = Year, y = Trend, group = Group)) +
    geom_line(alpha = 0.2, colour = "grey80") +
    geom_point(data = small, mapping = aes(x = Year, y = d15N), inherit.aes = FALSE,
               size = 0.7) +
    geom_line(data = newYear, mapping = aes(x = Year, y = fit), inherit.aes = FALSE) +
    labs(x = "Year", y = d15n_label)

plot_grid(plt1, plt2, plt3, ncol = 1, labels = "auto", align = "hv", axis = "lrtb",
          rel_widths = c(0.5, 1, 1))

