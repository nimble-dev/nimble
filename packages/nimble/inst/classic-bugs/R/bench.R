library(coda)
x <- read.openbugs(quiet=TRUE)
benchstats <- summary(x)$statistics
dump("benchstats", file=benchfile)
