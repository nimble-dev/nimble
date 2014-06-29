## By default we run examples in all sub-directories.  This can be
## overridden by the user-specified variable "EXAMPLES"

EXAMPLES ?= asia biops eyes stagnant pigs \
         alli beetles birats dugongs hearts jaw orange schools 

export JAGS ?= /usr/local/bin/jags

## Expand example names to make targets

bench_EX=$(addprefix bench_,$(EXAMPLES))
check_EX=$(addprefix check_,$(EXAMPLES))
Rcheck_EX=$(addprefix Rcheck_,$(EXAMPLES))
clean_EX=$(addprefix clean_,$(EXAMPLES))
distclean_EX=$(addprefix distclean_,$(EXAMPLES))

## All make targets are phony: i.e. they do not correspond to file
## names

.PHONY: all bench check Rcheck clean distclean \
	$(bench_EX) $(check_EX) $(Rcheck_EX) $(clean_EX) $(distclean_EX)

## Make targets

all:
bench: $(bench_EX)
check: $(check_EX)
Rcheck: $(Rcheck_EX)
clean: $(clean_EX)
distclean: $(distclean_EX)

## Dispatch targets to subdirectories. They can be run in parallel

$(bench_EX): 
	$(MAKE) -f $(CURDIR)/Makefile.sub -C $(patsubst bench_%,%,$@) bench 

$(check_EX): 
	$(MAKE) -f $(CURDIR)/Makefile.sub -C $(patsubst check_%,%,$@) check

$(Rcheck_EX): 
	$(MAKE) -f $(CURDIR)/Makefile.sub -C $(patsubst Rcheck_%,%,$@) Rcheck

$(clean_EX): 
	$(MAKE) -f $(CURDIR)/Makefile.sub -C $(patsubst clean_%,%,$@) clean

$(distclean_EX): 
	$(MAKE) -f $(CURDIR)/Makefile.sub -C $(patsubst distclean_%,%,$@) distclean
