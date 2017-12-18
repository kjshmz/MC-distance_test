TARGET = a.out

OBJECTS = main.o
MOD_FILES =
#FC = ifort -O2 -openmp -shared-intel -mcmodel=large #-parallel
#FC = /usr/local/openmpi-1.6.4_intel-13.1.3.192/bin/mpif90 -CB #-parallel #-shared-intel -mcmodel=large
FC = /usr/local/openmpi-1.5.1-intel64-v12.0.0u1/bin/mpif90 -shared-intel -mcmodel=large # -parallel
#FC = mpif90 -CB #-parallel #-shared-intel -mcmodel=large

FFLAGS =
LDFLAGS = #-openmp

#MKLROOT  = /opt/intel/composer_xe_2011/mkl
#FFLAGS  += -I${MKLROOT}/include/ia32 -I${MKLROOT}/include
#LDFLAGS += -L${MKLROOT}/lib/ia32 ${MKLROOT}/lib/ia32/libmklblas95.a
#LDFLAGS += ${MKLROOT}/lib/ia32/libmkl_lapack95.a
#LDFLAGS += -lmkl_intel -lmkl_intel_thread -lmkl_core -openmp -lpthread -lm

.SUFFIXES: .o .f90
.f90.o:
	${FC} -c $<

${TARGET}: ${OBJECTS}
	${FC} -o $@ ${OBJECTS} ${LDFLAGS} ${FFLAGS}

.PHONY: clean
clean:
	${RM} ${TARGET} ${OBJECTS} ${MOD_FILES}
