## Analyse the d15N data for Small Water, SMALL1 core

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
source("functions/simulate.gamm.R")

## Data
small <- readRDS("./data/small-water/small-water-isotope-data.rds")
head(small)

## Generate a plot of the data
plt <- ggplot(small, aes(x = Year, y = d15N)) +
    geom_point() +
        theme_bw() +
            ylab(expression(delta^{15}*N ~(Bulk~Organic~Matter)))
plt

## model it using gamm()
mod <- gamm(d15N ~ s(Year, k = 15), data = small,
            correlation = corCAR1(form = ~ Year), method = "REML")

summary(mod$gam)
intervals(mod$lme, which = "var-cov")
plot(mod$gam)

## Predict from model
newYear <- with(small, data.frame(Year = seq(min(Year), max(Year), length.out = 100)))
newYear <- cbind(newYear, data.frame(predict(mod$gam, newYear, se.fit = TRUE)))
crit <- qt(0.975, df = df.residual(mod$gam))
newYear <- transform(newYear,
                     upper = fit + (crit * se.fit),
                     lower = fit - (crit * se.fit))

## Draw the fitted spline on the data
plt.fit <- ggplot(small, aes(x = Year, y = d15N)) +
    geom_ribbon(data = newYear, aes(y = fit, ymax = upper, ymin = lower, x = Year),
                fill = "grey", alpha = 0.3) +
    geom_line(data = newYear, aes(y = fit, x = Year), col = "black", size = 0.8) +
    geom_point(col = "darkgrey") +
    theme_bw() +
    ylab(expression(delta^{15}*N ~(Bulk~Organic~Matter)))
plt.fit

## Need to do this manually, not through `fd` - write a function
set.seed(1)
nsim <- 50
randSims <- simulate(mod, nsim = nsim, newdata = newYear, unconditional = TRUE)
colnames(randSims) <- paste0("sim", seq_len(nsim))
randSims <- setNames(stack(as.data.frame(randSims)), c("simulated", "run"))
randSims <- transform(randSims, Year = rep(newYear$Year, nsim),
                      simulated = simulated)

plt.sim <- ggplot(newYear, aes(x = Year, y = fit)) +
    geom_line(data = randSims, mapping = aes(y = simulated, x = Year, group = run),
              col = "black", alpha = 0.1) +
    geom_line(size = 0.8, col = "black") +
    theme_bw() +
    ylab(expression(delta^{15}*N ~(Bulk~Organic~Matter))) +
    xlab("Year")
plt.sim

## Derivatives
fd <- derivSimulCI(mod, samples = 10000, n = 100)
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
    ylab(expression(italic(hat(f)*"'"*(Year)))) +
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

plt.fit2 <- plt.fit +
    geom_line(data = newYear, mapping = aes(y = YearIncr), size = 1.5) +
    geom_line(data = newYear, mapping = aes(y = YearDecr), size = 1.5)
plt.fit2

## Cowplot
displayPlt <- plot_grid(plt.fit, plt.sim, derivPlt, plt.fit2, labels = "AUTO")
displayPlt

ggsave("./figures/small-water-display.pdf", displayPlt)
