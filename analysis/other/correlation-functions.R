## Packages
library("ggplot2")
library("cowplot")
library("viridis")

## load the correlation functions
source("./functions/correlation-functions.R")

pExpdf <- expand.grid(h = seq(0, 7, length = 100),
                      phi = c(1,2),
                      k = c(0.5, 1, 2))
pExpdf <- transform(pExpdf, correlation = powerExp(h, phi, k),
                    kfac = factor(k))

pexp.plt <- ggplot(pExpdf, aes(x = h, y = correlation, group = kfac, colour = kfac)) +
    geom_line(size = 1) +
    facet_wrap( ~ phi, labeller = label_bquote(phi == .(phi))) +
    ylab(expression(rho)) +
    xlab(expression(italic(h))) +
    ##scale_color_viridis(option = "viridis", discrete = TRUE) +
    scale_colour_manual(values = viridis(11)[c(2,5,10)]) +
    theme_bw() +
    guides(colour = guide_legend(title = expression(kappa)))
pexp.plt

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
    ## scale_color_viridis(discrete = TRUE) +
    scale_colour_manual(values = viridis(11)[c(2,5,10)]) +
    theme_bw() +
    guides(colour = guide_legend(title = expression(kappa)))
pmat.plt

correlFuns <- plot_grid(pexp.plt, pmat.plt, ncol = 1, align = "hv",
                        labels = "AUTO")
correlFuns

ggsave("./figures/correlation-functions-plot.pdf", correlFuns, width = 9,
       height = 9)
