## Simple Age depth model for Small Water

## Should be run from base git repo

## Monotonic spline model

## R Packages
library("mgcv")                         # explicit but `monoSpline()` will check mgcv is available
library("ggplot2")                      # for plotting

source("./functions/monoSpline.R")

## SMALL1 core
small <- read.csv("./data/small-water/small1-dating.csv")
## monotonic spline age-depth model
smallAgeMod <- with(small, monoSpline(Depth, Date, Error))

## Predict from fitted monotonic spline
newDepths <- data.frame(Depth = seq(0, max(small$Depth), length.out = 100))
newDepths <- transform(newDepths, Fitted = predict(smallAgeMod, newdata = newDepths))

## Plot observed and interpolated ages
plt <- ggplot(small, aes(x = Date, y = Depth)) +
    geom_errorbarh(aes(xmin = Date - Error, xmax = Date + Error, height = 0)) +
        geom_point() +
            geom_line(data = newDepths, mapping = aes(x = Fitted, y = Depth), colour = "red") +
                scale_y_reverse() + scale_x_reverse() +
                    theme_bw() +
                        xlab("Year") + ylab("Sediment depth (cm)")
plt

## Add age info the Small Water Isotope Data
smallIso <- read.csv("./data/small-water/small-water-small1-isotope-data.csv")
smallIso <- transform(smallIso, Year = predict(smallAgeMod, newdata = smallIso))

## Save updated Small Water isotope data
saveRDS(smallIso, "./data/small-water/small-water-isotope-data.rds")
