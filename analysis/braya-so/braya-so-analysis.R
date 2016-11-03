## Analyse the UK37 data for Braya So core

## Packages
library("mgcv")
library("ggplot2")
library("cowplot")
library("viridis")

pal <- viridis(2)
theme_set(theme_bw())

## load the simultaneous confidence interval derivatives code
tmpf <- tempfile()
download.file("https://gist.githubusercontent.com/gavinsimpson/ca18c9c789ef5237dbc6/raw/295fc5cf7366c831ab166efaee42093a80622fa8/derivSimulCI.R",
              tmpf, method = "wget")
source(tmpf)

## Source other functions
source("functions/simulate.gam.R")

## Data
braya <- read.table("./data/braya-so/DAndrea.2011.Lake Braya So.txt",
                    skip = 84)
names(braya) <- c("Depth", "DepthUpper", "DepthLower", "Year", "YearYoung", "YearOld", "UK37")
braya <- transform(braya, sampleInterval = YearYoung - YearOld)
head(braya)

## Generate a plot of the data
ylabel <- expression(italic(U)[37]^{italic(k)})
plt <- ggplot(braya, aes(x = Year, y = UK37)) +
    geom_point(colour = pal[1]) +
    geom_line(colour = pal[1]) +
    theme_bw() +
    ylab(ylabel)
plt

## variography...
library("gstat")
library("sp")
m.reml <- gam(UK37 ~ s(Year, k = 30), data = braya, method = "REML")
resi.reml <- data.frame(Residuals = resid(m.reml), Year = braya$Year, Dummy = 1)
coordinates(resi.reml) <- ~ Year + Dummy
hscat(Residuals ~ 1, data = resi.reml, breaks = c(0, 25, 50, 100,200, 500, 1000))
resi.vgmCloud <- variogram(Residuals ~ 1, data = resi.reml, cloud = TRUE)
resi.vgm <- variogram(Residuals ~ 1, data = resi.reml, width = 20, cressie = TRUE)

ggplot(resi.vgm, aes(x = dist, y = gamma)) +
    geom_smooth(se = FALSE, col = "grey", span = 0.75,
                method.args = list(family = "gaussian")) +
    geom_point() +
    xlim(c(0, 1000)) + ylim(c(0, 1e-3)) +
    ylab(expression(Semivariance ~~ gamma)) +
    xlab("Temporal separation (years)")

## model it using gam()
effRange <- 250
mod <- gam(UK37 ~ s(Year, k = 30, bs = "gp", m = c(3, effRange)), data = braya,
           method = "REML")
summary(mod)
plot(mod, residuals = TRUE, pch = 16, n = 300)

mod2 <- gam(UK37 ~ s(Year, k = 30, bs = "gp", m = c(3, effRange)), data = braya,
            method = "REML", weights = sampleInterval)
summary(mod2)
plot(mod2, residuals = TRUE, pch = 16, n = 300)

mod3 <- gam(UK37 ~ s(Year, k = 30, bs = "ad"), data = braya,
            method = "REML", weights = sampleInterval)
summary(mod3)
plot(mod3, residuals = TRUE, pch = 16, n = 300)

## Model as a cr but spread knots at quantiles & use weights
knots <- with(braya, list(Year = unname(quantile(Year, probs = seq(0L, 1L, length.out = 35)))))
mod4 <- gam(UK37 ~ s(Year, k = 35, bs = "cr"), data = braya, knots = knots, method = "ML")
summary(mod4)
plot(mod4, residuals = TRUE, pch = 16, n = 300)
rug(knots[[1]], side = 3)

knots <- with(braya, list(Year = unname(quantile(Year, probs = seq(0L, 1L, length.out = 35)))))
mod5 <- gam(UK37 ~ s(Year, k = 35, bs = "cr"), data = braya, knots = knots, method = "REML",
            weights = sampleInterval, select = TRUE)
summary(mod5)
plot(mod5, residuals = TRUE, pch = 16, n = 300)
rug(knots[[1]], side = 3)

## Model 6: TPRS, weights as sampleInterval, k = needs to be higher
mod6 <- gam(UK37 ~ s(Year, k = 40, bs = "tp"), data = braya, method = "REML", weights = sampleInterval)
summary(mod6)
plot(mod6, residuals = TRUE, pch = 16, n = 300)

## Model 7
mod7 <- gam(UK37 ~ s(Year, k = 35, bs = "cr"), data = braya, method = "REML",
            weights = sampleInterval, select = TRUE)
summary(mod7)
plot(mod7, residuals = TRUE, pch = 16, n = 300)

## Model 8
mod8 <- gam(UK37 ~ s(Year, k = 40, bs = "ad"), data = braya, method = "REML")
summary(mod8)
plot(mod8, residuals = TRUE, pch = 16, n = 300)


## wrap this in a function that will return all the plots and derived objects
processGAM <- function(mod) {
    ## Predict from model
    N <- 500
    newYear <- with(braya, data.frame(Year = seq(min(Year), max(Year), length.out = N)))
    newYear <- cbind(newYear, data.frame(predict(mod, newYear, se.fit = TRUE)))

    ## Draw the fitted spline on the data
    plt.fit <- ggplot(braya, aes(x = Year, y = UK37)) +
        geom_line(data = newYear, aes(y = fit, x = Year), colour = pal[1]) +
            geom_point(colour = pal[1]) +
                theme_bw() +
                    ylab(ylabel)
    ## plt.fit

    ## Need to do this manually, not through `fd` - write a function
    set.seed(1)
    nsim <- 10000
    take <- 50
    sims <- simulate(mod, nsim = nsim, newdata = newYear, unconditional = TRUE)
    randSims <- sims[, sample(nsim, take)]
    colnames(randSims) <- paste0("sim", seq_len(take))
    randSims <- setNames(stack(as.data.frame(randSims)), c("simulated", "run"))
    randSims <- transform(randSims, Year = rep(newYear$Year, take),
                          simulated = simulated)

    ## Add a CI to the fitted spline
    CI <- apply(sims, 1L, quantile, probs = c(0.025, 0.975))
    newYear <- transform(newYear, lower = CI[1,], upper = CI[2,])

    plt.sim <- ggplot(newYear, aes(x = Year, y = fit)) +
        geom_line(data = randSims, mapping = aes(y = simulated, x = Year, group = run),
                  colour = pal[1], alpha = 0.1) +
            geom_line(colour = pal[1]) +
                theme_bw() +
                    ylab(ylabel) +
                        xlab("Year")
    ## plt.sim

    ## Derivatives
    fd <- derivSimulCI(mod, samples = 10000, n = N)
    CI <- apply(fd[[1]]$simulations, 1, quantile, probs = c(0.025, 0.975))
    sigD <- signifD(fd[["Year"]]$deriv, fd[["Year"]]$deriv, CI[2, ], CI[1, ],
                    eval = 0)
    newYear <- transform(newYear,
                         derivative = fd[["Year"]]$deriv[, 1], # computed first derivative
                         fdUpper = CI[2, ],                    # upper CI on first deriv
                         fdLower = CI[1, ],                    # lower CI on first deriv
                         increasing = sigD$incr,               # where is curve increasing?
                         decreasing = sigD$decr)               # ... or decreasing?

    ## Plot it
    derivPlt <- ggplot(newYear, aes(x = Year, y = derivative)) +
        geom_ribbon(aes(ymax = fdUpper, ymin = fdLower), alpha = 0.3, fill = pal[1]) +
            geom_line(colour = pal[1]) +
                geom_line(aes(y = increasing), colour = pal[1], size = 1.5) +
                    geom_line(aes(y = decreasing), colour = pal[1], size = 1.5) +
                        ylab(expression(italic(hat(f) * "'") * (Year))) +
                            xlab("Year") +
                                theme_bw()
    ## derivPlt

    ## Plot of data, fitted spline, CI, and significant periods of change
    ## Need the signifD helper function to generate the vectors for increasing
    ## and decreasing periods
    sigDYear <- signifD(newYear$fit, fd[["Year"]]$deriv[,1], CI[2, ], CI[1, ],
                        eval = 0)
    newYear <- transform(newYear,
                         YearIncr = sigDYear$incr,
                         YearDecr = sigDYear$decr)

    plt.fit2 <- ggplot(braya, aes(x = Year, y = UK37)) +
        geom_ribbon(data = newYear,
                    mapping = aes(ymax = upper, ymin = lower, y = fit, x = Year),
                    fill = pal[1], alpha = 0.3) +
            geom_line(data = newYear, aes(y = fit, x = Year), colour = pal[1]) +
                geom_point(colour = pal[1]) +
                    theme_bw() +
                        ylab(ylabel) +
                            geom_line(data = newYear, mapping = aes(y = YearIncr), size = 1.5, colour = pal[1]) +
                                geom_line(data = newYear, mapping = aes(y = YearDecr), size = 1.5, colour = pal[1])
    ## plt.fit2

    combinedPlt <- plot_grid(plt, plt.sim, derivPlt, plt.fit2, labels = "AUTO",
                             align = "vh")
    ## combinedPlt

    ## Return object:
    out <- list(plots = list(combinedPlt, plt.fit, plt.fit2, plt.sim, derivPlt),
                objects = newYear)
    out
}

pltsMod1 <- processGAM(mod = mod)       # Gaussian process smooth
pltsMod2 <- processGAM(mod = mod2)      # Gaussian process smooth with weights
pltsMod3 <- processGAM(mod = mod3)      # Adaptive smooth with weights
pltsMod4 <- processGAM(mod = mod4)      # Cubic spline with knots at quantiles
pltsMod5 <- processGAM(mod = mod5)      # Cubic spline with knots at quantiles & weights
pltsMod6 <- processGAM(mod = mod6)      # TP spline with weights
pltsMod7 <- processGAM(mod = mod7)      # Cubic spline with weights
pltsMod8 <- processGAM(mod = mod8)      # Adaptive spline

pltsMod1$plots[[1]]   # Gaussian process smooth
pltsMod2$plots[[1]]   # Gaussian process smooth with weights
pltsMod3$plots[[1]]   # Adaptive smooth with weights
pltsMod4$plots[[1]]   # Cubic spline with knots at quantiles
pltsMod5$plots[[1]]   # Cubic spline with knots at quantiles & weights
pltsMod6$plots[[1]]   # TP spline with weights
pltsMod7$plots[[1]]   # Cubic spline with weights
pltsMod8$plots[[1]]   # Adaptive spline

combinedPlt <- pltsMod2$plots[[1]]       # change this if needed to add sample weight information to GP smooth
ggsave("./figures/braya-so-display.pdf", combinedPlt, width = 10, height = 7)

pltData <- do.call("rbind", lapply(list(pltsMod1, pltsMod2, pltsMod3, pltsMod6, pltsMod7, pltsMod8), `[[`, "objects"))
pltData <- transform(pltData, Model = rep(c("GP", "GP+weights", "AD+weights", "TPRS+weights", "CRS+weights", "AD"),
                              each = nrow(pltsMod1$objects)))

allFits <- ggplot(pltData, aes(x = Year, y = fit)) +
    geom_point(aes(x = Year, y = UK37), data = braya) +
    geom_line(aes(colour = Model)) + labs(y = ylabel, x = "Year")
allFits

ggsave("./figures/braya-so-fits-using-various-spline-bases.pdf", allFits, height = 4, width = 9)
ggsave("./figures/braya-so-fits-using-various-spline-bases.png", allFits)
