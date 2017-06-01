# this Makefile needs to be cleaned up and tested - CJP 6/29/14

nimble_0.1.tar.gz :  nimble/*
	R CMD build nimble

prep:  nimble/R nimble_0.1.tar.gz nimble/NAMESPACE nimble/man
      # need installed pkg to be able to list all fxns/classes
      R CMD INSTALL nimble_0.1.tar.gz
      ./prep_pkg.R

all:  nimble/man/* nimble_0.1.tar.gz nimble/NAMESPACE

clean :
	rm -f nimble_0.1.tar.gz
	# add rm -f nimble/man* eventually
