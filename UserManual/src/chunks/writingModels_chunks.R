## @knitr linearRegressionGraph
library(nimble)
mc <- nimbleCode({
    intercept ~ dnorm(0, sd = 1000)
    slope ~ dnorm(0, sd = 1000)
    sigma ~ dunif(0, 100)
    for(i in 1:4) {
        predicted.y[i] <- intercept + slope * x[i]
        y[i] ~ dnorm(predicted.y[i], sd = sigma)
    }
})

model <- nimbleModel(mc, data = list(y = rnorm(4)))

library(igraph)

layout <- matrix(ncol = 2, byrow = TRUE,
   # These seem to be rescaled to fit in the plot area,
   # so I'll just use 0-100 as the scale
                 data = c(33, 100,
                          66, 100,
                          50, 0, # first three are parameters
                          15, 50, 35, 50, 55, 50, 75, 50, # x's
                          20, 75, 40, 75, 60, 75, 80, 75, # predicted.y's
                          25, 25, 45, 25, 65, 25, 85, 25) # y's
                 )

sizes <- c(45, 30, 30,
           rep(20, 4),
           rep(50, 4),
           rep(20, 4))

edge.color <- "black"
    # c(
    # rep("green", 8),
    # rep("red", 4),
    # rep("blue", 4),
    # rep("purple", 4))
stoch.color <- "deepskyblue2"
det.color <- "orchid3"
rhs.color <- "gray73"
fill.color <- c(
    rep(stoch.color, 3),
    rep(rhs.color, 4),
    rep(det.color, 4),
    rep(stoch.color, 4)
)


plot(model$graph, vertex.shape = "crectangle",
     vertex.size = sizes,
     vertex.size2 = 20,
     layout = layout,
     vertex.label.cex = 3.0,
     vertex.color = fill.color,
     edge.width = 3,
     asp = 0.5,
     edge.color = edge.color)

## @knitr linearRegressionCode
{
    intercept ~ dnorm(0, sd = 1000)
    slope ~ dnorm(0, sd = 1000)
    sigma ~ dunif(0, 100)
    for(i in 1:4) {
        predicted.y[i] <- intercept + slope * x[i]
        y[i] ~ dnorm(predicted.y[i], sd = sigma)
    }
}

## @knitr linearRegressionAltCode
{
    intercept ~ dnorm(0, sd = 1000)
    slope ~ dnorm(0, sd = 1000)
    sigma ~ dunif(0, 100)
    for(i in 1:4) {
        y[i] ~ dnorm(intercept + slope * x[i], sd = sigma)
    }
}


