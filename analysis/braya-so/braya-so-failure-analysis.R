## Analyse the UK37 data for Braya So core
##
## Illustrates how f(x) + CAR(1) + e model fails

## Packages
library("mgcv")
library("ggplot2")
library("cowplot")
library("viridis")

## load the simultaneous confidence interval derivatives code
tmpf <- tempfile()
download.file("https://gist.githubusercontent.com/gavinsimpson/ca18c9c789ef5237dbc6/raw/295fc5cf7366c831ab166efaee42093a80622fa8/derivSimulCI.R",
              tmpf, method = "wget")
source(tmpf)

## Source other functions
source("functions/simulate.gamm.R")

## Data
braya <- read.table("./data/braya-so/DAndrea.2011.Lake Braya So.txt",
                    skip = 84)
names(braya) <- c("Depth", "DepthUpper", "DepthLower", "Year", "YearYoung", "YearOld", "UK37")
head(braya)

mod <- gamm(UK37 ~ s(Year, k = 20), data = braya,
            correlation = corCAR1(form = ~ Year),
            method = "REML")
summary(mod$gam)
intervals(mod$lme)
plot(mod$gam, resid = TRUE, pch = 16)

## ggplot with data and fitted spline, then resids vs time in second panel
ylabel <- expression(italic(U)[37]^{italic(k)})
N <- 300
newYear <- with(braya, data.frame(Year = seq(min(Year), max(Year), length.out = N)))
newYear <- cbind(newYear, data.frame(predict(mod$gam, newYear, se.fit = TRUE)))
crit.t <- qt(0.975, df = df.residual(mod$gam))
newYear <- transform(newYear,
                     upper = fit + (crit.t * se.fit),
                     lower = fit - (crit.t * se.fit))

mod2 <- gam(UK37 ~ s(Year, k = 30), data = braya)
summary(mod2)
plot(mod2, n = N)
gam.check(mod2)

fit <- predict(mod2, newdata = newYear, se.fit = TRUE)

newYear <- rbind(newYear, newYear)
newYear[seq(N+1, length.out = N, by = 1), ]$fit <- fit$fit
newYear[seq(N+1, length.out = N, by = 1), ]$upper <-
    fit$fit + (qt(0.975, df.residual(mod2)) * fit$se.fit)
newYear[seq(301, length.out = 300, by = 1), ]$lower <-
    fit$fit - (qt(0.975, df.residual(mod2)) * fit$se.fit)
newYear <- transform(newYear, Method = rep(c("GAMM (REML)", "GAM (GCV)"), each = N))

p1 <- ggplot(braya, aes(y = UK37, x = Year)) +
    geom_point() +
    geom_ribbon(data = newYear,
                mapping = aes(y = fit, x = Year, ymax = upper, ymin = lower, fill = Method),
                alpha = 0.4) +
    geom_line(data = newYear, mapping = aes(y = fit, x = Year, colour = Method), size = 1) +
    ylab(ylabel) +
    theme_bw() +
    scale_color_manual(values = viridis(11)[c(2, 6)]) +
    scale_fill_manual(values = viridis(11)[c(2, 6)]) +
    theme(legend.position = "top")
p1

tmp <- cbind(rbind(braya, braya),
             residuals = c(resid(mod$lme), resid(mod2, type = "response")),
             Method = rep(c("GAMM (REML)", "GAM (GCV)"), each = nrow(braya)))
p2 <- ggplot(tmp,
             aes(y = residuals, x = Year, colour = Method)) +
    geom_point() +
    geom_line() +
    ylab("Normalised residuals") +
    theme_bw() +
    scale_color_manual(values = viridis(11)[c(2, 6)]) +
    theme(legend.position = "top")
p2

combined <- plot_grid(p1, p2, labels = "AUTO", align = "vh")
combined

ggsave("./figures/braya-so-failure-figure.pdf", width = 9, height = 4)
