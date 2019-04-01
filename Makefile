FC = f95
#FCFLAGS = -O0 -g -fbounds-check
FCFLAGS = -O0 -g 
PROGRAMS = debug
RLIBS = tree_climb.so

all: $(PROGRAMS) $(RLIBS)

debug: tree_climb.o

%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)
%.o: %.f90
	$(FC) $(FCFLAGS) -c $<
.PHONY: clean
clean:
	rm -f *.o *.so *~ $(PROGRAMS)
%.so: %.f90
	R CMD SHLIB tree_climb.f90
