## Analyse the UK37 data for Braya So core

## Packages
library("mgcv")
library("ggplot2")
library("cowplot")

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
head(braya)

## Generate a plot of the data
ylabel <- expression(italic(U)[37]^{italic(k)})
plt <- ggplot(braya, aes(x = Year, y = UK37)) +
    geom_point(col = "darkgrey") +
    geom_line(col = "darkgrey") +
    theme_bw() +
    ylab(ylabel)
plt

## model it using gam()
mod <- gam(UK37 ~ s(Year, k = 50, bs = "gp", m = c(3, 100)), data = braya, method = "REML")
summary(mod)
plot(mod)

## Predict from model
N <- 300
newYear <- with(braya, data.frame(Year = seq(min(Year), max(Year), length.out = N)))
newYear <- cbind(newYear, data.frame(predict(mod, newYear, se.fit = TRUE)))

## Draw the fitted spline on the data
plt.fit <- ggplot(braya, aes(x = Year, y = UK37)) +
    geom_line(data = newYear, aes(y = fit, x = Year), col = "black") +
    geom_point(col = "darkgrey") +
    theme_bw() +
    ylab(ylabel)
plt.fit

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
              col = "black", alpha = 0.1) +
    geom_line() +
    theme_bw() +
    ylab(ylabel) +
    xlab("Year")
plt.sim

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
    geom_ribbon(aes(ymax = fdUpper, ymin = fdLower), alpha = 0.3, fill = "grey") +
    geom_line() +
    geom_line(aes(y = increasing), col = "black", size = 1.5) +
    geom_line(aes(y = decreasing), col = "black", size = 1.5) +
    ylab(expression(italic(hat(f) * "'") * (Year))) +
    xlab("Year") +
    theme_bw()
derivPlt

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
                fill = "grey", alpha = 0.3) +
    geom_line(data = newYear, aes(y = fit, x = Year), col = "black") +
    geom_point(col = "darkgrey") +
    theme_bw() +
    ylab(ylabel) +
    geom_line(data = newYear, mapping = aes(y = YearIncr), size = 1.5) +
    geom_line(data = newYear, mapping = aes(y = YearDecr), size = 1.5)
plt.fit2

combinedPlt <- plot_grid(plt, plt.sim, derivPlt, plt.fit2, labels = "AUTO",
                         align = "vh")
combinedPlt

ggsave("./figures/braya-so-display.pdf", combinedPlt, width = 10, height = 7)
