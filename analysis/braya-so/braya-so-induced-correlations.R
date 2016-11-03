# ---- induced-correlations-plt ----

## Analyse the UK37 data for Braya So core: Look at induced correlations

## Packages
library("mgcv")
library("ggplot2")
library("cowplot")
library("grid") # for unit.pmax(), unit.list()

theme_set(theme_bw())

## Source other functions
source("functions/induced-correlations.R")

## Data
braya <- read.table("./data/braya-so/DAndrea.2011.Lake Braya So.txt",
                    skip = 84)
names(braya) <- c("Depth", "DepthUpper", "DepthLower", "Year", "YearYoung", "YearOld", "UK37")
braya <- transform(braya, sampleInterval = YearYoung - YearOld)

## induced correlations....
T <- 601
t <- 1:T
K <- 10
set.seed(1)
et <- rnorm(T, sd = 0.25)
ft <- sin(2 * pi * ((t-1) / T))
yt <- ft + et
egData <- data.frame(yt = yt, t = t, ft = ft)

egData.plt <- ggplot(egData, aes(x = t, y = yt)) +
    geom_point(size = 0.9) +
        geom_line(aes(y = ft), col = "red", lwd = 1) +
            labs(x = expression(italic(t)), y = expression(y[italic(t)]))
## egData.plt

m.ps <- gam(yt ~ s(t, bs = "ps", k = K), data = egData, method = "REML")
m.tp <- gam(yt ~ s(t, bs = "tp", k = K), data = egData, method = "REML")

tp <- c(100,225,350,500)

egDataCor <- rbind(addInduced(t, tp, m.ps), addInduced(t, tp, m.tp))
egDataCor <- cbind(egDataCor, Basis = rep(c("P Spline", "Thinplate Spline"), each = T * length(tp)))

eg.plt <- ggplot(egDataCor, aes(x = t, y = correl, colour = tp)) +
    geom_line() + facet_wrap( ~ Basis) + scale_colour_discrete(guide = FALSE) +
            labs(x = expression(italic(t)), y = "Induced correlation")
## eg.plt

uk.ps <- gam(UK37 ~ s(Year, k = 30, bs = "ps"), data = braya, method = "REML")
uk.tp <- gam(UK37 ~ s(Year, k = 30, bs = "tp"), data = braya, method = "REML")

ukTP <- c(20,40,60,80)
ukData <- with(braya, rbind(addInduced(Year, ukTP, uk.ps), addInduced(Year, ukTP, uk.tp)))
ukData <- cbind(ukData, Basis = rep(c("P Spline", "Thinplate Spline"), each = nrow(braya) * length(ukTP)))

uk.plt <- ggplot(ukData, aes(x = t, y = correl, colour = tp)) +
    geom_line() + facet_wrap( ~ Basis) +
        labs(x = "Year BP", y = "Induced correlation") + scale_colour_discrete(guide = FALSE)
## uk.plt

## cowplot the panels together
## Need to sort out the widths of the plot elements on the left iof each plot
g.egData <- ggplotGrob(egData.plt) # convert to gtable
g.eg <- ggplotGrob(eg.plt)         # convert to gtable
g.uk <- ggplotGrob(uk.plt)         # convert to gtable
egData.widths <- g.egData$widths[1:3]   # extract the first three widths,
                                        # corresponding to left margin, y lab, and y axis
eg.widths <- g.eg$widths[1:3] # same for mpg plot
uk.widths <- g.uk$widths[1:3] # same for mpg plot
max.widths <- unit.pmax(egData.widths, eg.widths, uk.widths) # calculate maximum widths
g.egData$widths[1:3] <- max.widths  # assign max widths to egData gtable
g.eg$widths[1:3] <- max.widths      # assign max widths to eg.plt gtable
g.uk$widths[1:3] <- max.widths      # assign max widths to uk.plt gtable

induced.plt <- plot_grid(g.egData, g.eg, g.uk, ncol = 1, rel_heights = c(1, 1, 1), labels = "AUTO")

## ggsave("./figures/induced-correlations.pdf", plot = induced.plt, height = 8, width = 8)

## plot it
induced.plt
