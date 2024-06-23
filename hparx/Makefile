# Makefile to make HPARX library

INC_COMMON = ./Common.mk
include $(INC_COMMON)

## Macros
TARGETS = libhparx.a
LIB_SRCS = globals.f90 hparx_base.f90 hparx_math.f90 hparx_mtab.f90 \
	hparx_file.f90 hparx_rand.f90 hparx_radi.f90 hparx_rayl.f90 \
	hparx_geo3D.f90 hparx_lina.f90 hparx_stat.f90 hparx_vecmat.f90 \
	hparx.f90 hparx_nons.f90
#	hparx_iMie.f90 hparx_sscat.f90 hparx_atmos.f90 \

LIB_OBJS = $(LIB_SRCS:.f90=.o)
BUPIDX  = $(shell date "+%Y%m%d_%H%M")
BUPFILE = hparx-src.$(BUPIDX).tar.gz

## Commands & rules
all : $(TARGETS)

clean :
	rm -f *.o *.mod $(TARGETS) *~
echo :
	@ echo $(LIB_SRCS)
print :
	@ cat $(LIB_SRCS)
backup : clean
	@ rm -f $(BUPFILE)
	tar cvzf $(BUPFILE) *.f90 Makefile *.mk
.f90.o : 
	$(FC) -c $< $(FCFLAGS)

## Targets
libhparx.a : $(LIB_OBJS)
	@ rm -f $@
	@ echo "Archiving..." $@
	@ if test $(RANLIB); then \
		$(AR) cr $@ $^; ranlib $@; \
        else \
		$(AR) crs $@ $^; \
	fi
	@ echo "Done."

## Dependency
hparx.o : globals.o hparx_base.o hparx_math.o hparx_mtab.o hparx_vecmat.o \
	hparx_file.o hparx_lina.o hparx_stat.o hparx_rand.o hparx_radi.o \
	hparx_rayl.o hparx_geo3D.o hparx_nons.o
#hparx_iMie.o hparx_sscat.o hparx_atmos.o \


#hparx_atmos.o  : globals.o hparx_math.o
hparx_base.o   : globals.o
hparx_file.o   : globals.o hparx_base.o
hparx_geo3D.o  : globals.o hparx_math.o hparx_rand.o
#hparx_iMie.o   : globals.o hparx_math.o
hparx_lina.o   : globals.o hparx_base.o hparx_vecmat.o
hparx_math.o   : globals.o hparx_lina.o hparx_vecmat.o
hparx_mtab.o   : globals.o
hparx_rand.o   : globals.o hparx_base.o
hparx_radi.o   : globals.o hparx_math.o hparx_rand.o hparx_stat.o hparx_vecmat.o
hparx_rayl.o   : globals.o hparx_rand.o
#hparx_sscat.o  : globals.o hparx_base.o hparx_math.o hparx_file.o
hparx_stat.o   : globals.o hparx_math.o hparx_lina.o
hparx_vecmat.o : globals.o
hparx_nons.o   : globals.o hparx_base.o
