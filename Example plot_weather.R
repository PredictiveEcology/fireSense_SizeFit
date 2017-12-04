library(dplyr)
library(ggplot2)
library(PtProcess)
library(SpaDES)

set.seed(1)

modulePath <- "~/Documents/GitHub/McIntire-lab/modulesPrivate/"
start <- end <- 1

# Define simulation parameters
times <- list(start = start, end = end, timeunit = "year")
modules <- list("fireSense_SizeFit", "fireSense_SizePredict")
paths <- list(
  modulePath = modulePath
)

# Model fitting: Create random weather and fire size dataset
dummyData <- data.frame(
  weather = scale(rep(1:100, each = 10)),
  fireSize = rtappareto(1000, lambda = rep(seq(3, .3, length.out = 100), each = 10), theta = rep(seq(100, 2e4, length.out = 100), each = 10), a = 1)
)

# Predict: Create a data.frame with the covariate(s)
predData <- data.frame( weather = (seq(1, 100, length.out = 10) - mean(1:100)) / sd(rep(1:100, each = 10)))


# Define module parameters
parameters <- list(
  fireSense_SizeFit = list(
    formula = list(beta = fireSize ~ weather,
                   theta = fireSize ~ weather),
    a = 1,
    data = "dummyData"
  ),
  fireSense_SizePredict = list(
    modelName = "fireSense_SizeFitted",
    data = "predData"
  )
)

# Objects to pass from the global environment to the simList environment
objects <- c("dummyData", "predData")

# Create the simList
sim <- simInit(
  times = times, 
  params = parameters, 
  modules = modules, 
  objects = objects, 
  paths = paths
)

sim <- spades(sim)

# Plots
# Influence of weather of the fire size distribution
sizeSeq <- seq(1, 1e5, length.out = 1000)
probs <- with(sim$fireSense_SizePredicted[[as.character(start)]],
              unlist(Map(f = ptappareto, lambda = beta, theta = theta, MoreArgs = list(q = sizeSeq, a = a, lower.tail = FALSE))))

plotData <- data_frame(prob = pmax(probs, .Machine$double.neg.eps), group = rep(1:10, each = length(sizeSeq)))

x11()
cols <- colorRampPalette(c("#FED976", "#FEB24C", "#FD8D3C", "#FC4E2A", "#E31A1C", "#BD0026", "#800026"))(10)
p <- ggplot(data = plotData) + theme_bw()
p <- p + coord_cartesian(xlim = range(sizeSeq), ylim = c(1e-4, 1), expand = FALSE)
p <- p + geom_line(aes(x = rep(sizeSeq, 10), y = prob, group = group, color = group), size = 2)
p <- p + scale_x_log10(name = "Fire size (ha)", breaks = 10^(0:5), labels = format(10^(0:5), big.mark = " ", trim = TRUE, scientific = FALSE)) 
p <- p + scale_y_log10(name = "Probability of reaching a size > X ha", breaks = c(1e-4, 1e-3, 1e-2, .1, 1), labels = c("0.01%", "0.1%", "1%", "10%", "100%")) 
p <- p + scale_colour_gradientn(space = "Lab", colours = cols, breaks = c(1, 10), limits = c(1, 10), expand = c(0,0),
                                guide = guide_colorbar(direction = "horizontal", title = "Drought intensity", title.position = "top", label = TRUE,
                                                       barwidth = unit(6, "lines"), barheight = unit(1.6, "lines"), draw.llim = TRUE, draw.ulim = TRUE))
p


# Predictions against observations
# Transform raw x,y coordinates to correct x,y coordinates for plotting survival functions
survX <- function(x) c(x[1], rep(x[-1], each = 2))
survY <- function(y) c(rep(y[-length(y)], each = 2), y[length(y)])

pp <- ppoints(sort(dummyData$fireSize))

surv_tapPareto <- bind_rows(
  !!! lapply(1:100, function(i) with(sim$fireSense_SizePredicted[[as.character(start)]],
                                     quantile(rtappareto(1000, lambda = rep(beta, each = 100), theta = rep(theta, each = 100), a = a), 1-pp)) )
)

plot(survX(quantile(dummyData$fireSize, 1-pp)), survY(pp), type = "l", log = "xy", xlab = "Fire size", ylab = "Probability of reaching a size > X ha")
for (i in 1:100) lines(survX(unlist(surv_tapPareto[i,])), survY(pp), col = gray(.8))
lines(survX(quantile(dummyData$fireSize, 1-pp)), survY(pp), col = "red", lwd = 2)


