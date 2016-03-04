## Analyse the d15N data for Small Water, SMALL1 core

## Packages
library("mgcv")
library("ggplot2")

## load the simultaneous confidence interval derivatives code
tmpf <- tempfile()
download.file("https://gist.githubusercontent.com/gavinsimpson/ca18c9c789ef5237dbc6/raw/295fc5cf7366c831ab166efaee42093a80622fa8/derivSimulCI.R",
              tmpf, method = "wget")
source(tmpf)

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
mod <- gamm(d15N ~ s(Year, k = 15), data = small, correlation = corCAR1(form = ~ Year), method = "REML")

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
            geom_line(data = newYear, aes(y = fit, x = Year), col = "darkgrey", size = 0.8) +
                geom_point() +
        theme_bw() +
            ylab(expression(delta^{15}*N ~(Bulk~Organic~Matter)))
plt.fit

## Derivatives
fd <- derivSimulCI(mod, samples = 10000, n = 100)
CI <- apply(fd[[1]]$simulations, 1, quantile, probs = c(0.025, 0.975))

## Need to do this manually, not through `fd` - write a function
set.seed(1)
nsim <- 20
randSims <- fd[[1]]$simulations[, sample(10000, nsim)]
colnames(randSims) <- paste0("sim", 1:20)
randSims <- setNames(stack(as.data.frame(randSims)), c("simulated", "run"))
randSims <- transform(randSims, Year = rep(newYear$Year, nsim),
                      simulated = simulated + coef(mod$gam)[1])

plt.sim <- ggplot(newYear, aes(x = Year, y = fit)) +
    geom_line(size = 1) +
        geom_line(data = randSims, mapping = aes(y = simulated, x = Year, group = run), col = "grey") +
        theme_bw()
plt.sim

