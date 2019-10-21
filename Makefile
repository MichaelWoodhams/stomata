# When installing the R library, builds will be done using R library
# default build rules. This Makefile is just to build the 'debug'
# test program, suitable for singlestepping with a debugger.
# It is quick and simple rather than carefully crafted to allow for
# future expansion.

FC = f95
FCFLAGS = -O0 -g -fbounds-check
PROGRAM = debug_tree_climb
SRCDIR = fasttree2019/src
LIB = $(SRCDIR)/fasttree2019.so
LIBSOURCE = $(SRCDIR)/fasttree2019.f95

all: $(LIB) $(PROGRAM) 

$(PROGRAM): $(PROGRAM).f95 $(LIB)
	$(FC) $(FCFLAGS) $^ -o $(PROGRAM)

clean:
	rm -f $(SRCDIR)/*.o $(SRCDIR)/*.so $(PROGRAM)

$(LIB): $(LIBSOURCE)
	R CMD SHLIB $(LIBSOURCE)
