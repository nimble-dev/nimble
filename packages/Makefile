.PHONY: FORCE

VERSION=$(shell grep Version: nimble/DESCRIPTION | grep -o '[^ ]*$$')

version: FORCE
	@echo $(VERSION)

configure: nimble/configure.ac
	(cd nimble; autoconf)

# this rebuilds Rd files in man/ and rebuilds NAMESPACE
man: FORCE
	./prep_pkg.R

nimble_$(VERSION).tar.gz: nimble/configure FORCE
	@echo "Building nimble package in this directory"
	(R CMD build nimble)

build: nimble_$(VERSION).tar.gz

install: nimble_$(VERSION).tar.gz
	R CMD INSTALL nimble_$(VERSION).tar.gz

check: build FORCE
	R CMD check --as-cran nimble_$(VERSION).tar.gz

rhub: build FORCE
	Rscript -e "library(rhub); check(platform = c('ubuntu-gcc-devel','windows-x86_64-devel','linux-x86_64-rocker-gcc-san','linux-x86_64-centos6-epel','fedora-clang-devel'))"

test: FORCE
	cd .. ; ./install_requirements.R
	cd .. ; ./run_tests.R --parallel

benchmark: FORCE
	mkdir -p profile
	R --quiet --slave --vanilla --file=benchmark.R

profile: FORCE
	mkdir -p profile
	NIMBLE_BENCHMARK_DIR=$(PWD)/profile \
	NIMBLE_BENCHMARK_SEC=10.0 \
	  R --quiet --slave --vanilla --file=benchmark.R \
	   --debugger=valgrind \
	   --debugger-args="--tool=callgrind --callgrind-out-file=profile/callgrind.out"
	kcachegrind profile/callgrind.out || callgrind_annotate profile/callgrind.out

# We are progressively delegating style to clang-format.
# Some older C++ files are still styled manually.
clang-format: FORCE
	clang-format -style=file \
	  nimble/inst/include/nimble/NimArr.h \
	  nimble/inst/include/nimble/NimArrBase.h \
	  nimble/inst/include/nimble/predefinedNimbleLists.h \
	  nimble/inst/include/nimble/tensorflow.h \
	  nimble/inst/CppCode/predefinedNimbleLists.cpp \
	  nimble/inst/CppCode/tensorflow.cpp \
	  -i

# Keeps .Rproj.user/ and .Rhistory.
clean: FORCE
	git clean -dfx -e .R*

mrproper: FORCE
	git clean -dfx

FORCE:

