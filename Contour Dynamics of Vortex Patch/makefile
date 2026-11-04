# 
#	Fortran 90 rules
#
.SUFFIXES: .f90 .m $(SUFFIXES)
#
FC90=f95
LIB=$(FORT_LIB)
# F90FLAGS=-O3  -qalias=noaryovrlp
COMPILE.f90=$(FC90) $(F90FLAGS) -c
LINK.f90=$(FC90) $(LNKFLAGS) $(F90FLAGS)

.f90.o: 
	$(COMPILE.f90) $<
.f90.m:
	$(COMPILE.f90) $<
.f90:
	$(FC90) $(F90FLAGS) $(LNKFLAGS) -o $@ $< $(LIBS)



# dependencies

dump_.m:      dump_.f90 dyn_com.m
cubic_.m:     cubic_.f90 dyn_com.m
dmarray_.m:   dmarray_.f90 dyn_com.m
uv_.m:        uv_.f90 cubic_.m 
advect_.m:    advect_.f90 uv_.m 
surgery_.m:   surgery_.f90 dmarray_.m dyn_com.m 
renode_.m:    renode_.f90 cubic_.m dmarray_.m
dyn.o:          dyn.f90 advect_.m dump_.m renode_.m \
		tabgen_.m	
#
SOURCES=advect_.m cubic_.m dmarray_.m surgery_.m\
    dump_.m dyn.o dyn_com.m renode_.m \
    uv_.m tabgen_.m 
#
dyn: $(SOURCES)
	$(LINK.f90) -o $@ $(SOURCES:.m=.o)

