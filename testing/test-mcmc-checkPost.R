# write code to simulate data from the model and see if posterior contains the true parameters
dataNodes <- Cmodel$getNodeNames(dataOnly = TRUE)

# need to set reasonable init values for each model, then simulate from data nodes


consts <- list(
           'pump' = list(N = 10, 
             t = c(94.3, 15.7, 62.9, 126, 5.24, 31.4, 1.05, 1.05, 2.1, 10.5)
             ),
           'blocker' = list( )
           )

trueParams <- list(
                'pump' = list(alpha = 1, beta = 1)
                )

example = 'pump'; vol = 'vol1'

dir <- ""
if(!"package:nimble" %in% search()) {
  dir <- file.path("..", "..", "packages", "nimble", "inst", "classic-bugs")
  dir <- file.path(dir, vol, example)
}

Rmodel <- readBUGSmodel(example, dir = dir, inits = trueParams[[example]], data = consts[[example]], useData = FALSE, useInits = FALSE)
simulate(Rmodel)

for(node in Rmodel$getNodeNames(dataOnly = TRUE))
  Rmodel$setData(node)
# now need to do setData; how determine bottom stoch nodes?


