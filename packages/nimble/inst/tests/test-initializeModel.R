source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))

RwarnLevel <- options('warn')$warn
options(warn = 1)
nimbleVerboseSetting <- nimbleOptions('verbose')
nimbleOptions(verbose = FALSE)

context("Testing of model initializaion using initializeModel")


test_that('initializeModel works', {
    code <- nimbleCode({
        x[1] <- 0.5
        x[2] <- x[1] + 1
        x[3] <- x[2]^2
        x[4] ~ dnorm(x[3], 1)
        x[5] <- x[4] + 1
        y ~ dnorm(x[5], 1)
        x[7] <- y + 1
        x[8] <- x[7] + 10
        x[9] ~ dnorm(x[8], 1)
        x[10] ~ dnorm(x[9], 1)
    })

    Rmodel <- nimbleModel(code, data=list(y = 5), calculate = FALSE)

    Rinit <- initializeModel(Rmodel)
    Cmodel <- compileNimble(Rmodel)
    Cinit <- compileNimble(Rinit, project = Rmodel)

    correctValues <-  c(0.5, 1.5, 2.25, 3.512954, 4.512954, NA, 6, 16, 15.673767, 17.003566)

    set.seed(0)
    Rinit$run()
    expect_identical(as.numeric(round(Rmodel$x, 6)), correctValues)
    expect_true(is.na(Rmodel$x[6]))

    set.seed(0)
    Cinit$run()
    expect_identical(as.numeric(round(Cmodel$x, 6)), correctValues)
    expect_true(is.na(Cmodel$x[6]))
})




test_that('initializeModel works correctly for state space models', {
    set.seed(0)
    nCues=2
    nTrials=44*nCues
    validity=0.75
    t=1
    changeT<-round(runif(nCues-1,40,48))
    changeT<-cumsum(changeT)
    genTarget=sign(runif(nTrials,-10000,10000))
    genTarget[genTarget==0]<-1
    relCue=rep(NA,nTrials)
    opt=c(1,2,3)
    relCue[1]=sample(opt,size = 1)
    genCues<-matrix(rep(NA,nTrials*3),nrow=nTrials)
    genCues[1,relCue[1]]<-genTarget[1]
    other=c(1,2,3)
    other<-other[other!=relCue[1]]
    genCues[1,other]<-sample(c(-1,1),1)
    genCues[1,genCues[1,]==0]<-1
    ##
    for (t in 2:nTrials){
        relCue[t]=relCue[t-1]
        if (max(t==changeT)==1) {
            if (length(opt)>1) {relCue[t]=sample(opt,size = 1)} else {relCue[t]=opt}
        }
        valid=runif(1,0,1)
        if (valid<=validity){genCues[t,relCue[t]]<-genTarget[t]} else {genCues[t,relCue[t]]<--1*genTarget[t]}
        genCues[t,-relCue[t]]<-sign(runif(2,-10000,10000))
    }
    ##
    test=data.frame(genCues,genTarget,relCue)
    ##
    stateSpaceCode<-nimbleCode({
        gamma<-0.75
        tau<-0.9
        x[1]~dcat(lambda[1,1:3])
        lambda[1,1]<-1/3
        lambda[1,2]<-1/3
        lambda[1,3]<-1/3
        prob[1]<-equals(cues[1,x[1]],1)*gamma+equals(cues[1,x[1]],-1)*(1-gamma)
        target[1]~dbern(prob[1])
        for (t in 2:nTrials){
            x[t]~dcat(lambda[t,1:3])
            lambda[t,1]<-equals(x[t-1],1)*tau+(1-equals(x[t-1],1))*(1-tau)/2
            lambda[t,2]<-equals(x[t-1],2)*tau+(1-equals(x[t-1],2))*(1-tau)/2
            lambda[t,3]<-equals(x[t-1],3)*tau+(1-equals(x[t-1],3))*(1-tau)/2
            prob[t]<-equals(cues[t,x[t]],1)*gamma+equals(cues[t,x[t]],-1)*(1-gamma)
            target[t]~dbern(prob[t])
        }
    })
    ##
    datalist<-list(cues=as.matrix(test[,1:3]),target=as.numeric((test[,4]==1)))
    inits<-list(gamma=0.75, tau=0.9)
    constants<-list(nTrials=length(datalist$target))
    
    Rmodel<-nimbleModel(code = stateSpaceCode,data = datalist,constants=constants,inits = inits, calculate = FALSE)
    
    Rinit <- initializeModel(Rmodel)
    Cmodel <- compileNimble(Rmodel)
    Cinit <- compileNimble(Rinit, project = Rmodel)
    
    set.seed(0)
    Rinit$run()
    
    set.seed(0)
    Cinit$run()
    
    expect_equal(Rmodel$calculate(), -77.49817, tol = 0.000001)
    expect_equal(Cmodel$calculate(), -77.49817, tol = 0.000001)
    
    stateSpaceModel<-nimbleModel(code = stateSpaceCode,data = datalist,constants=constants,inits = inits, calculate = FALSE)
    bootstrapFilter<-buildBootstrapFilter(stateSpaceModel,nodes='prob',control=list(saveAll=TRUE))
    compiledList<-compileNimble(stateSpaceModel,bootstrapFilter)
    
    set.seed(0)
    expect_equal(bootstrapFilter$run(10), -52.78133, tol = 0.000001)
    
    set.seed(0)
    expect_equal(compiledList$bootstrapFilter$run(10), -52.78133, tol = 0.000001)
})



test_that('initializeModel recalculates all deterministic nodes in topological order', {
    set.seed(0)
    code <- nimbleCode({
        p <- ilogit(a)
        x ~ dbern(p)
    })
    constants <- list()
    data <- list()
    inits <- list(a = 10)
    Rmodel <- nimbleModel(code, constants, data, inits)
    Cmodel <- compileNimble(Rmodel)
    ##
    expect_true(is.na(Rmodel$calculate()))
    expect_true(is.na(Cmodel$calculate()))
    ##
    my_initializeModel <- initializeModel(Rmodel)
    Cmy_initializeModel <- compileNimble(my_initializeModel, project = Rmodel)
    ##
    my_initializeModel$run()
    Cmy_initializeModel$run()
    ##
    expect_equal(Rmodel$x, 1)
    expect_equal(Cmodel$x, 1)
    ##
    expect_equal(Rmodel$getLogProb(), -0.0000453989, tol = 0.0000001)
    expect_equal(Cmodel$getLogProb(), -0.0000453989, tol = 0.0000001)
    ##
    Rmodel$x <- NA
    Cmodel$x <- NA
    Rmodel$a <- -10
    Cmodel$a <- -10
    ##
    my_initializeModel$run()
    Cmy_initializeModel$run()
    ##
    expect_equal(Rmodel$x, 0)
    expect_equal(Cmodel$x, 0)
    expect_equal(Rmodel$getLogProb(), -0.0000453989, tol = 0.0000001)
    expect_equal(Cmodel$getLogProb(), -0.0000453989, tol = 0.0000001)
})



options(warn = RwarnLevel)
nimbleOptions(verbose = nimbleVerboseSetting)


