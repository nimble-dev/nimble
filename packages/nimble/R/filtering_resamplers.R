resamplerVirtual <- nimbleFunctionVirtual(
  run = function(wts = integer(1)) {
    returnType(double(1))
  }
)

multinomialResampleFunction <- nimbleFunction(
  setup = function(){},
  run = function(wts = double(1)){
    n <- size(wts)
    idCounts <- rmulti(1, n, wts)
    idCounter <- 1
    for(i in 1:n){
      for(j in 1:idCounts[i]){
        ids[idCounter] <- i
      }
    }
    returnType(double(1))
    return(ids)
  },
  contains = resamplerVirtual
)


residualResampleFunction <- nimbleFunction(
  setup = function(){},
  run = function(wts = double(1)){
    n <- length(wts)
    if(n == 1){
      return(numeric(1, 1))
    }
    ids <- integer(n, 0)
    sumIds <- numeric(n, 0)
    expectedN <- wts*n
    floorN <- floor(expectedN)
    sumFloorN <- sum(floorN)
    remainder <- n - sumFloorN
    if(remainder > 0){
      barWts <- (expectedN - floorN)/remainder
      sumIds <- rmulti(1, remainder,  barWts)
      numberOfSamples <- floorN + sumIds
      idCounter <- 1
      for(i in 1:n){
        if(numberOfSamples[i] > 0){
          ids[idCounter:(idCounter + numberOfSamples[i] - 1)] <- i
          idCounter <- idCounter + numberOfSamples[i]
        }
      }
    }
    else{
      idCounter <- 1
      for(i in 1:n){
        if(floorN[i] > 0){
          ids[idCounter:(idCounter + floorN[i] - 1)] <- i
          idCounter <- idCounter + floorN[i]
        }
      }
    }
    returnType(integer(1))
    return(ids)
  },
  contains = resamplerVirtual
)

stratifiedResampleFunction <- nimbleFunction(
  setup = function(){},
  run = function(wts = double(1)){
    n <- length(wts)
    if(n == 1){
      return(numeric(1, 1))
    }
    ids <- integer(n)
    randUs <- numeric(n)
    sumWts <- numeric(n+1)
    for(i in 1:n){
      randUs[i] <- (runif(1, 0, 1) + (i - 1))/n
      sumWts[i+1] <- sumWts[i] + wts[i]
    }
    wtsCounter <- 2
    for(i in 1:n){
      while(randUs[i] >= sumWts[wtsCounter]){
        wtsCounter <- wtsCounter + 1
      }
      ids[i] <- wtsCounter - 1
    }
    returnType(integer(1))
    return(ids)
  },
  contains = resamplerVirtual
)

systematicResampleFunction <- nimbleFunction(
  setup = function(){},
  run = function(wts = double(1)){
    n <- length(wts)
    if(n == 1){
      return(numeric(1, 1))
    }
    ids <- integer(n)
    sumWts <- numeric(n)
    sumWts[1] <- wts[1]
    for(i in 2:n){
      sumWts[i] <- sumWts[i-1] + wts[i]
    }
    unifRV <- runif(1, 0, 1/n)
    wtsCounter <- 1
    for(i in 1:n){
      while(unifRV >= sumWts[wtsCounter]){
        wtsCounter <- wtsCounter + 1
      }
      ids[i] <- wtsCounter
      unifRV <- unifRV + 1/n
    }
    returnType(integer(1))
    return(ids)
  },
  contains = resamplerVirtual
)
