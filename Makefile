PROG=nmaker

HEADERS=dbginfo.h domain.h dtime.h global.h initial.h iosnap.h nmk_util.h offload.h partition.h partmesh.h ppkernel.h proto.h queue.h stepping.h subtree.h traversal.h
OBJS_CPU=domain.o global.o initial.o iosnap.o main.o partition.o partmesh.o ppkernel.o queue.o stepping.o subtree.o traversal.o
OBJS_MIC=domain.om global.om initial.om iosnap.om main.om partition.om partmesh.om ppkernel.om queue.om stepping.om subtree.om traversal.om

#OPENMP=-openmp
OPENMP=
#PAPI_PATH=/opt/experf/papi-5.3.2
#PAPI=-D__PAPI -I$(PAPI_PATH)/include $(PAPI_PATH)/lib/libpapi.a
#PAPI=
KMP_AFFINITY ?= compact
NUM_PROCS ?= 4
#NO_OFFLOAD=
NO_OFFLOAD=-no-offload
STATIC=
#STATIC=-static_mpi

MKLROOT?=/opt/intel/mkl

#PRECISION=-DNMK_DOUBLE_PREC
#FFTDEF=-DFFTW3_LIB
DEBUG=-DNMK_DEBUG
NAIVE=-DNMK_NAIVE_GRAVITY
#VERIFY=-DNMK_VERIFY #-DNMK_GEN_VER_FILE

CC=mpiicc
MICCC=icc
CFLAGS=$(STATIC) -xHost -fimf-domain-exclusion=15 -O3 $(DEBUG) $(NAIVE) $(VERIFY) $(PRECISION) $(OPENMP) $(NO_OFFLOAD) $(FFTDEF) -g -traceback -Bdynamic -fno-omit-frame-pointer 
CFLAGSMIC=-mmic -fimf-domain-exclusion=15 -O3 $(DEBUG) $(NAIVE) $(VERIFY) $(PRECISION) $(OPENMP) -g -traceback -Bdynamic -fno-omit-frame-pointer

INC=-I$(MKLROOT)/include/fftw

#FFTLIB=-L$(MKLROOT)/lib/intel64 -lfftw3x_cdft_lp64 -lfftw3xc_double_intel -lmkl_cdft_core -lmkl_blacs_intelmpi_lp64 -lmkl_core -lmkl_intel_lp64 -lmkl_sequential 
FFTLIB=-L$(MKLROOT)/lib/intel64 -lfftw2x_cdft_DOUBLE_lp64 -lmkl_cdft_core -lmkl_blacs_intelmpi_lp64 -lmkl_core -lmkl_intel_lp64 -lmkl_sequential 

LD=mpiicc
LDFLAGS=$(STATIC) $(FFTLIB) -traceback -Bdynamic -fno-omit-frame-pointer


all: $(PROG)

run: $(PROG)
	OMP_NUM_THREADS=$(OMP_NUM_THREADS) \
	KMP_AFFINITY=$(KMP_AFFINITY) \
	mpirun -np $(NUM_PROCS) ./$(PROG)

$(PROG): $(OBJS_CPU)
	$(LD) -o $(PROG) $(OBJS_CPU) $(LDFLAGS) $(OPENMP) $(NO_OFFLOAD)

pp:
	$(CC) -openmp -o pp9 pp9.c ppkernel.c $(CFLAGS) $(PAPI) $(NO_OFFLOAD)

ppm:
	$(CC) -openmp -o pp9m pp9m.c ppkernel.c $(CFLAGS) $(PAPI) $(NO_OFFLOAD)

ppmic: ppkernel.om
	$(MICCC) -openmp -o pp9.mic pp9.c ppkernel.om $(CFLAGSMIC)

ppmmic: ppkernel.om
	$(MICCC) -openmp -o pp9m.mic pp9m.c ppkernel.om $(CFLAGSMIC)

ppv:
	$(CC) -o ppverify ppverify.c 

.PHONY:clean
clean:
	rm -rf *.o *.om $(PROG) pp9 pp9.mic pp9m pp9m.mic *.dbg

.SUFFIXES: .c .o .om

.c.o:
	$(CC) $(CFLAGS) $(INC) -c $<

.c.om:
	$(MICCC) $(CFLAGSMIC) $(INC) -openmp -o $*.om -c $<

domain.o: domain.c global.h offload.h proto.h nmk_utils.h dbginfo.h dtime.h domain.h
global.o: global.c global.h offload.h proto.h nmk_utils.h dbginfo.h dtime.h
initial.o: initial.c global.h offload.h proto.h nmk_utils.h dbginfo.h dtime.h initial.h iosnap.h
iosnap.o: iosnap.c global.h offload.h proto.h nmk_utils.h dbginfo.h dtime.h iosnap.h
main.o: main.c global.h offload.h proto.h nmk_utils.h dbginfo.h dtime.h initial.h partition.h domain.h partmesh.h subtree.h traversal.h queue.h stepping.h
partition.o: partition.c global.h offload.h proto.h nmk_utils.h dbginfo.h dtime.h partition.h
partmesh.o: partmesh.c global.h offload.h proto.h nmk_utils.h dbginfo.h dtime.h partmesh.h
ppkernel.o: ppkernel.c ppkernel.h global.h offload.h proto.h nmk_utils.h dbginfo.h dtime.h
queue.o: queue.c global.h offload.h proto.h nmk_utils.h dbginfo.h dtime.h queue.h
stepping.o: stepping.c global.h offload.h proto.h nmk_utils.h dbginfo.h dtime.h stepping.h
subtree.o: subtree.c global.h offload.h proto.h nmk_utils.h dbginfo.h dtime.h domain.h subtree.h
traversal.o: traversal.c global.h offload.h proto.h nmk_utils.h dbginfo.h dtime.h domain.h queue.h ppkernel.h traversal.h

ppkernel.om: ppkernel.c ppkernel.h offload.h proto.h nmk_utils.h dbginfo.h dtime.h
