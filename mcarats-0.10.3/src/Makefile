# User macros
## MPI? {1 if using MPI}
MPI_USE = 0
## Suffix for MPI binary executables
MPI_SFX =

# Include common macros
INC_COMMON = ./Common.mk
include ${INC_COMMON}

# System macros
BUPIDX  = ${shell date "+%Y%m%d_%H%M"}
BUPFILE = mcarats-0.10.3-src.${BUPIDX}.tar.gz

# Files for mcarats
TARGETS = mcarats
MCAR_SRCS1 = mcarUtl.f90 mcarDsm.f90 mcarRpv.f90 mcarLsrt.f90 mcarPho.f90 mcarVis.f90 \
	mcarSrc.f90 main_mcarats.f90
MCAR_SRCS2 = mcarAtm.F90 mcarSfc.F90 mcarSca.F90 mcarFlx.F90 mcarRad.F90 mcarWld.F90 \
	mcarats.F90 globals_MPI.F90
MCAR_OBJS = ${MCAR_SRCS1:.f90=.o} ${MCAR_SRCS2:.F90=.o} ${LIB_OBJS}

# Commands
all : ${TARGETS}

install : all
	@ mkdir -p ${BINDIR}
	cp ${TARGETS} ${BINDIR}
clean :
	rm -f *.o *.mod ${TARGETS} *~
backup : clean
	pushd ..
	@ rm -f ${BUPFILE}
	tar -Lcvzf ${BUPFILE} ${MCAR_SRC1} ${MCAR_SRC2} Makefile *.mk COPYING
	popd
help :
	@ echo "usage: make [TARGET] [MACRO=value ...]"
	@ echo " "
	@ echo "TARGET"
	@ echo "  all     : compile all"
	@ echo "  install : install all"
	@ echo "  clean   : remove object files & executables"
	@ echo "  backup  : make an archive file of all important files"
	@ echo "  help    : print usage"
	@ echo "MACRO"
	@ echo "  {name=default : explanations}"
	@ echo "  MPI_USE=0     : 1 if using MPI"
	@ echo "  MPI_SFX=      : suffix for MPI executables"


# Suffix rules
.F90.o :
	${FC} -c $< ${FCFLAGS} ${LIB_FLAGS} ${FC_DEF}UseMPI=${MPI_USE}
.f90.o :
	${FC} -c $< ${FCFLAGS} ${LIB_FLAGS} 

# Targets
mcarats : ${MCAR_OBJS}
	@ echo Linking...
	@ if test ${MPI_USE} = 1; then \
		${FC} -o mcarats${MPI_SFX} $^ ${FCFLAGS} ${FC_MPI}; \
	else \
		${FC} -o mcarats${MPI_SFX} $^ ${FCFLAGS}; \
	fi
	@ echo Done: $@

# Dependency
main_mcarats.o : ${LIB_OBJS} mcarats.o
mcarats.o  : ${LIB_OBJS} globals_MPI.o mcarWld.o
mcarWld.o  : ${LIB_OBJS} globals_MPI.o mcarUtl.o mcarPho.o mcarSca.o mcarAtm.o \
             mcarSfc.o mcarSrc.o mcarFlx.o mcarRad.o mcarVis.o
mcarFlx.o  : ${LIB_OBJS} globals_MPI.o mcarUtl.o mcarPho.o
mcarRad.o  : ${LIB_OBJS} globals_MPI.o mcarUtl.o mcarPho.o mcarAtm.o mcarSfc.o \
             mcarVis.o
mcarAtm.o  : ${LIB_OBJS} globals_MPI.o mcarUtl.o mcarPho.o mcarSca.o
mcarSfc.o  : ${LIB_OBJS} globals_MPI.o mcarUtl.o mcarDsm.o mcarRpv.o mcarLsrt.o
mcarSrc.o  : ${LIB_OBJS} mcarAtm.o mcarSfc.o
mcarVis.o  : ${LIB_OBJS} mcarAtm.o mcarSfc.o
mcarSca.o  : ${LIB_OBJS} mcarUtl.o
mcarPho.o  : ${LIB_OBJS}
mcarDsm.o  : ${LIB_OBJS} mcarUtl.o
mcarRpv.o  : ${LIB_OBJS} mcarUtl.o
mcarLsrt.o : ${LIB_OBJS} mcarUtl.o
mcarUtl.o  : ${LIB_OBJS}
