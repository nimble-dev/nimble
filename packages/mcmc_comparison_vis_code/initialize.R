allModels <- c(# vol1
  'blocker',
  'bones',
  'dyes',
  'equiv', 
  #'line', 
  #'oxford', 
  'pump',
  'rats'
  # vol2
  #,'dugongs'
)

dir='~/Documents/R Stuff/Comparison'

# load all functions
sapply(list.files(pattern="*.R", path=paste(dir,"/functions",sep=''), full.names=TRUE), source);

# run mcmc suite for all examples, outputs are saved in the list 'mods'. each list element contains parameters and time
mods=runComparison(allModels,mcmcs=c('nimble','jags','noConj' ,'stan'),niter=10000,thin=1,
                   standir='~/Documents/example-models/bugs_examples/vol1/')

# makes html pages using mcmc outputs
# home page is called 'main.html'
make_pages(mods,allModels,ncol=3,dir)
